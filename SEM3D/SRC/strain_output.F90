!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file strain_output.F90
!!\brief Sortie des deformations a chaque pas de temps, pour le calcul des noyaux
!!       de sensibilite (inversion par etat adjoint).
!!
!! Active par "save_strain = true;" dans input.spec. Seuls les elements solides hors PML
!! (DM_SOLID_CG) sont ecrits, element par element (sans assemblage : la deformation est
!! discontinue entre elements). Un rang MPI ecrit deux fichiers dans ./strain :
!!  - strain_geom.RRRR.h5 : coords (npts,3), weight (npts) = w_i w_j w_k |J|, elem (npts),
!!    attributs dt, ngll, npts ;
!!  - strain.RRRR.bin : flux binaire float64, un enregistrement (npts,6) par pas de temps,
!!    composantes xx, yy, zz, xy, xz, yz du tenseur des deformations.
!! Comme pour les traces, l''enregistrement k (k = 0, 1, ...) contient le champ au temps (k+1)*dt.
!! Avec "strain_period = p;", seul un pas sur p est ecrit : l''enregistrement j contient le champ
!! au temps (j+1)*p*dt (attribut period du fichier de geometrie).
!!
!! "strain_mode = m;" choisit ce qui est ecrit dans ./strain (meme geometrie) :
!!  - m = 0 : les deformations, comme ci-dessus ;
!!  - m = 1 (simulation directe) : disp.RRRR.bin, le deplacement aux noeuds des elements ecrits
!!    (3 composantes par noeud, liste des noeuds fixe), un enregistrement par periode ; 6 fois
!!    moins de donnees que les deformations, qui en sont recalculees a l''identique ;
!!  - m = 2 (simulation adjointe) : lit ./fwd_strain/disp.RRRR.bin, ecrit par la simulation
!!    directe (meme maillage, meme decoupage), recalcule les deformations directes et accumule
!!    les noyaux en chaque point (enregistrement adjoint r <-> direct M-2-r, M enregistrements) :
!!      K_lambda = - p dt sum tr eps(u) tr eps(w),  K_mu = - 2 p dt sum eps(u) : eps(w),
!!    ecrits a la fin dans kernel.RRRR.bin (K_lambda (npts) puis K_mu (npts), float64).
!!    Les deformations adjointes ne sont pas ecrites.
!!
!! "strain_pml = true;" (modes 1 et 2) ajoute les points des PML solides (DM_SOLID_CG_PML), apres
!! ceux du domaine solide, dans la geometrie et les noyaux (attribut npts_pml). Dans la PML, la
!! deformation est la deformation etiree eps~ = C^-1 : (sigma1 + sigma2 + sigma3), qui donne la
!! contrainte totale (somme des contraintes scindees, amorties) : les noyaux y suivent les memes
!! formules avec eps~. Au bord interieur de la PML (amortissement nul), eps~ = eps(u).
!! En frequence, le noyau PML est S eps~(u) : dC : eps~(w), S = sx sy sz, sk = 1 + dk/(i omega),
!! dk l''amortissement dans la direction k : en temps, S f = f + a1 I1 + a2 I2 + a3 I3, Ik l''integrale
!! k-ieme de f en temps, a1 = dx+dy+dz, a2 = dx dy+dx dz+dy dz, a3 = dx dy dz. La simulation directe
!! ecrit S eps~(u) aux points PML dans epml.RRRR.bin (float32, (npts_pml,6) par enregistrement,
!! integrales sur les enregistrements), que la simulation adjointe relit dans ./fwd_strain.
!<
#include "index.h"
#include "optims.h"

module mstrain_output
    use sdomain
    use constants
    use semdatafiles
    use sem_c_bindings, only : sem_mkdir
    implicit none
    private
    public :: init_strain_output, save_strain_output, close_strain_output

    integer, save :: strain_unit
    logical, save :: strain_open = .false.
    integer, save :: npts = 0
    integer, save :: period = 1, nstep = 0
    integer, allocatable, dimension(:), save :: elems, pelems
    integer, save :: npts_s = 0
    real(kind=8), allocatable, dimension(:,:), save :: buffer
    integer, save :: mode = 0, nnod = 0, nrecf = 0, disp_unit, rank_of_kernels = 0
    integer, allocatable, dimension(:), save :: nodes      ! noeuds du domaine solide ecrits
    integer, allocatable, dimension(:), save :: nodemap    ! indice du domaine -> rang dans nodes
    real(kind=8), allocatable, dimension(:,:), save :: dispf, bufferf, ufull
    real(kind=8), allocatable, dimension(:), save :: kl, km
    integer, save :: epml_unit
    real(kind=4), allocatable, dimension(:,:), save :: epf     ! eps~ aux points PML (float32)
    real(kind=8), allocatable, dimension(:,:), save :: acoef   ! a1, a2, a3 aux points PML
    real(kind=8), allocatable, dimension(:,:,:), save :: iint  ! I1, I2, I3 de eps~ (directe)
    ! bloc et indice dans le bloc (lnum/VCHUNK, mod(lnum,VCHUNK)) des elements ecrits, solides et PML :
    ! les boucles de calcul (sur GPU avec OpenACC) n''accedent pas a Tdomain%specel
    integer, allocatable, dimension(:), save :: elb, ele, pelb, pele

contains

    subroutine strain_filenames(rank, fgeom, fdata)
        integer, intent(in) :: rank
        character(len=MAX_FILE_SIZE), intent(out) :: fgeom, fdata
        fgeom = "strain/strain_geom."//trim(adjustl(strrank(rank)))//".h5"
        fdata = "strain/strain."//trim(adjustl(strrank(rank)))//".bin"
    end subroutine strain_filenames

    subroutine init_strain_output(Tdomain)
        use HDF5
        use sem_hdf5
        use sem_mpi
        type(domain), intent(inout) :: Tdomain
        !
        integer :: n, nel, i, j, k, ip, ngll, bnum, ee, ierr, hdferr
        integer(HID_T) :: fid
        integer(kind=8) :: fsize
        character(len=MAX_FILE_SIZE) :: fgeom, fdata
        real(fpp), allocatable, dimension(:,:) :: coords
        real(fpp), allocatable, dimension(:) :: weight
        integer, allocatable, dimension(:) :: elem_of_pt
        real(fpp), dimension(:), allocatable :: tmp

        if (Tdomain%rank == 0) then
            ierr = sem_mkdir("strain")
        end if
        call MPI_Barrier(Tdomain%communicateur, ierr)

        mode = Tdomain%config%strain_mode
        call list_elems(Tdomain, DM_SOLID_CG, elems)
#ifdef CPML
        if (Tdomain%config%strain_pml /= 0) stop "strain_pml : PML classiques seulement (OPT_CPML = OFF)"
#endif
        if (Tdomain%config%strain_pml /= 0 .and. mode /= 0) then
            call list_elems(Tdomain, DM_SOLID_CG_PML, pelems)
        else
            allocate(pelems(0:-1))
        end if
        nel = size(elems)
        call block_of(Tdomain, elems, elb, ele)
        call block_of(Tdomain, pelems, pelb, pele)

        ngll = Tdomain%sdom%ngll
        npts_s = nel*ngll*ngll*ngll
        npts = npts_s + size(pelems)*ngll*ngll*ngll
        allocate(buffer(0:5, 0:max(npts,1)-1))
        allocate(coords(0:2, 0:max(npts,1)-1), weight(0:max(npts,1)-1), elem_of_pt(0:max(npts,1)-1))

        ip = 0
        call geom_of(Tdomain, Tdomain%sdom, elems, ip, coords, weight, elem_of_pt)
        call geom_of(Tdomain, Tdomain%spmldom, pelems, ip, coords, weight, elem_of_pt)

        call strain_filenames(Tdomain%rank, fgeom, fdata)
        call h5fcreate_f(trim(fgeom), H5F_ACC_TRUNC_F, fid, hdferr)
        if (npts > 0) then
            call write_dataset_d2(fid, "coords", coords)
            call write_dataset_d1(fid, "weight", weight)
            allocate(tmp(0:npts-1))
            tmp = real(elem_of_pt, fpp)
            call write_dataset_d1(fid, "elem", tmp)
            deallocate(tmp)
        end if
        call write_attr_real(fid, "dt", Tdomain%TimeD%dtmin)
        call write_attr_int(fid, "ngll", ngll)
        call write_attr_int(fid, "npts", npts)
        call write_attr_int(fid, "npts_pml", npts - npts_s)
        period = max(1, Tdomain%config%strain_period)
        nstep = 0
        call write_attr_int(fid, "period", period)
        call h5fclose_f(fid, hdferr)
        deallocate(coords, weight, elem_of_pt)

        rank_of_kernels = Tdomain%rank
        if (mode == 0) then
            open(newunit=strain_unit, file=trim(fdata), access="stream", form="unformatted", &
                status="replace", action="write")
        else
            call init_nodes(Tdomain)
            if (mode == 1) then
                open(newunit=disp_unit, file="strain/disp."//trim(adjustl(strrank(Tdomain%rank)))//".bin", &
                    access="stream", form="unformatted", status="replace", action="write")
            else
                call init_kernels(Tdomain)
            end if
            if (npts > npts_s) then
                allocate(epf(0:5, 0:npts-npts_s-1))
                if (mode == 1) then
                    allocate(acoef(0:2, 0:npts-npts_s-1), iint(0:5, 0:npts-npts_s-1, 1:3))
                    iint = 0d0
                    call pml_coefs(Tdomain)
                end if
                if (mode == 1) then
                    open(newunit=epml_unit, file="strain/epml."//trim(adjustl(strrank(Tdomain%rank)))//".bin", &
                        access="stream", form="unformatted", status="replace", action="write")
                else
                    open(newunit=epml_unit, file="fwd_strain/epml."//trim(adjustl(strrank(Tdomain%rank)))//".bin", &
                        access="stream", form="unformatted", status="old", action="read")
                    if (nnod == 0) then  ! rang sans element solide : enregistrements comptes sur epml
                        inquire(unit=epml_unit, size=fsize)
                        nrecf = int(fsize / (int(npts-npts_s, 8) * 24_8))
                    end if
                end if
            end if
        end if
        ! Donnees sur le GPU pour toute la simulation ; seuls le deplacement aux noeuds ecrits (dispf),
        ! eps~ aux points PML (epf) et, a la fin, les noyaux (kl, km) passent par l''hote.
        !$acc enter data copyin(elb, ele, pelb, pele) create(buffer)
        if (mode /= 0) then
            !$acc enter data copyin(nodes) create(dispf)
            if (npts > npts_s) then
                !$acc enter data create(epf)
                if (mode == 1) then
                    !$acc enter data copyin(acoef, iint)
                end if
            end if
            if (mode == 2) then
                !$acc enter data copyin(kl, km, ufull) create(bufferf)
            end if
        end if
        strain_open = .true.
        if (Tdomain%rank == 0) write(*,*) "--> STRAIN OUTPUT: ", nel, " + ", size(pelems), &
            " PML elements, ", npts, " points on rank 0"
    end subroutine init_strain_output

    !> Elements du domaine dm, dans l''ordre de Tdomain%specel.
    subroutine list_elems(Tdomain, dm, el)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: dm
        integer, allocatable, dimension(:), intent(out) :: el
        integer :: n, nel
        nel = count(Tdomain%specel(0:Tdomain%n_elem-1)%domain == dm)
        allocate(el(0:nel-1))
        nel = 0
        do n = 0, Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain == dm) then
                el(nel) = n
                nel = nel + 1
            end if
        end do
    end subroutine list_elems

    !> Bloc (b) et indice dans le bloc (e) de chaque element de el.
    subroutine block_of(Tdomain, el, b, e)
        type(domain), intent(in) :: Tdomain
        integer, dimension(0:), intent(in) :: el
        integer, allocatable, dimension(:), intent(out) :: b, e
        integer :: n
        allocate(b(0:max(size(el),1)-1), e(0:max(size(el),1)-1))
        b = 0
        e = 0
        do n = 0, size(el)-1
            b(n) = Tdomain%specel(el(n))%lnum/VCHUNK
            e(n) = mod(Tdomain%specel(el(n))%lnum, VCHUNK)
        end do
    end subroutine block_of

    !> Coordonnees, poids et element des points des elements el du domaine dom, a partir de ip.
    subroutine geom_of(Tdomain, dom, el, ip, coords, weight, elem_of_pt)
        type(domain), intent(in) :: Tdomain
        class(dombase), intent(in) :: dom
        integer, dimension(0:), intent(in) :: el
        integer, intent(inout) :: ip
        real(fpp), dimension(0:,0:), intent(inout) :: coords
        real(fpp), dimension(0:), intent(inout) :: weight
        integer, dimension(0:), intent(inout) :: elem_of_pt
        integer :: n, i, j, k, ngll, bnum, ee
        ngll = dom%ngll
        do n = 0, size(el)-1
            bnum = Tdomain%specel(el(n))%lnum/VCHUNK
            ee = mod(Tdomain%specel(el(n))%lnum, VCHUNK)
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        coords(:,ip) = Tdomain%GlobCoord(:, Tdomain%specel(el(n))%Iglobnum(i,j,k))
                        weight(ip) = dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(i,j,k,bnum,ee)
                        elem_of_pt(ip) = el(n)
                        ip = ip + 1
                    end do
                end do
            end do
        end do
    end subroutine geom_of

    !> Noeuds (indices du domaine dom) des elements el, et l''application inverse map.
    subroutine nodes_of(Tdomain, dom, el, nd, map)
        type(domain), intent(in) :: Tdomain
        class(dombase), intent(in) :: dom
        integer, dimension(0:), intent(in) :: el
        integer, allocatable, dimension(:), intent(out) :: nd, map
        integer :: n, i, j, k, ngll, bnum, ee, ind, nn
        ngll = dom%ngll
        allocate(map(0:max(dom%nglltot,0)))
        map = -1
        nn = 0
        do n = 0, size(el)-1
            bnum = Tdomain%specel(el(n))%lnum/VCHUNK
            ee = mod(Tdomain%specel(el(n))%lnum, VCHUNK)
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        if (map(ind) < 0) then
                            map(ind) = nn
                            nn = nn + 1
                        end if
                    end do
                end do
            end do
        end do
        allocate(nd(0:nn-1))
        do ind = 0, ubound(map, 1)
            if (map(ind) >= 0) nd(map(ind)) = ind
        end do
    end subroutine nodes_of

    subroutine init_nodes(Tdomain)
        type(domain), intent(inout) :: Tdomain
        call nodes_of(Tdomain, Tdomain%sdom, elems, nodes, nodemap)
        nnod = size(nodes)
        allocate(dispf(0:2, 0:max(nnod,1)-1))
    end subroutine init_nodes

    subroutine init_kernels(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer(kind=8) :: fsize
        character(len=MAX_FILE_SIZE) :: fname
        fname = "fwd_strain/disp."//trim(adjustl(strrank(Tdomain%rank)))//".bin"
        inquire(file=trim(fname), size=fsize)
        if (fsize < 0) stop "strain_mode = 2 : fwd_strain/disp.RRRR.bin introuvable"
        nrecf = 0
        if (nnod > 0) nrecf = int(fsize / (int(nnod, 8) * 24_8))
        open(newunit=disp_unit, file=trim(fname), access="stream", form="unformatted", &
            status="old", action="read")
        allocate(kl(0:max(npts,1)-1), km(0:max(npts,1)-1))
        kl = 0d0
        km = 0d0
        allocate(bufferf(0:5, 0:max(npts,1)-1))
        allocate(ufull(0:Tdomain%sdom%nglltot, 0:2))
        ufull = 0d0
    end subroutine init_kernels

    !> Deformations aux points du domaine solide (out(:,0:npts_s-1)) du champ de deplacement U
    !! (indice du domaine solide). U et out sont sur le GPU avec OpenACC.
    subroutine strain_of(Tdomain, U, out)
        type(domain), intent(inout) :: Tdomain
        real(fpp), dimension(0:,0:), intent(in) :: U
        real(kind=8), dimension(0:,0:), intent(inout) :: out
        call strain_sol(Tdomain%sdom, size(elems), U, out)
    end subroutine strain_of

    !> Deformation etiree eps~ = C^-1 : sigma aux points PML (out(:,npts_s:npts-1)), sigma etant la
    !! somme des contraintes scindees courantes (xx, yy, zz, xy, xz, yz).
    subroutine strain_pml_of(Tdomain, out)
        type(domain), intent(inout) :: Tdomain
        real(kind=8), dimension(0:,0:), intent(inout) :: out
#ifndef CPML
        call strain_spml(Tdomain%spmldom, size(pelems), npts_s, out)
#endif
    end subroutine strain_pml_of

#ifndef CPML
    subroutine strain_spml(dom, npel, off, out)
        type(domain_solidpml), intent(in) :: dom
        integer, intent(in) :: npel, off
        real(kind=8), dimension(0:,0:), intent(inout) :: out
        integer :: n, i, j, k, ngll, bnum, ee, ip, c
        real(kind=8) :: sig(0:5), lambda, mu, tr
        ngll = dom%ngll
        !$acc parallel loop collapse(4) async(1) firstprivate(ngll, off) private(bnum, ee, ip, c, sig, lambda, mu, tr) &
        !$acc& present(dom, dom%m_Stress1, dom%m_Stress2, dom%m_Stress3, dom%m_Lambda, dom%m_Mu, pelb, pele, out)
        do n = 0, npel-1
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        bnum = pelb(n)
                        ee = pele(n)
                        ip = off + ((n*ngll + k)*ngll + j)*ngll + i
                        do c = 0, 5
                            sig(c) = dom%m_Stress1(ee,c,i,j,k,bnum) + dom%m_Stress2(ee,c,i,j,k,bnum) &
                                + dom%m_Stress3(ee,c,i,j,k,bnum)
                        end do
                        lambda = dom%Lambda_(i,j,k,bnum,ee)
                        mu = dom%Mu_(i,j,k,bnum,ee)
                        tr = (sig(0) + sig(1) + sig(2)) / (3d0*lambda + 2d0*mu)
                        do c = 0, 2
                            out(c,ip) = (sig(c) - lambda*tr) / (2d0*mu)
                        end do
                        do c = 3, 5
                            out(c,ip) = sig(c) / (2d0*mu)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine strain_spml
#endif

    !> a1, a2, a3 aux points PML, a partir des amortissements dk = 2 (1/DumpSk(1) - 1) / dt.
    subroutine pml_coefs(Tdomain)
        type(domain), intent(inout) :: Tdomain
#ifndef CPML
        integer :: n, i, j, k, ngll, bnum, ee, ip
        real(kind=8) :: dx, dy, dz, c
        ngll = Tdomain%spmldom%ngll
        c = 2d0 / Tdomain%TimeD%dtmin
        ip = 0
        do n = 0, size(pelems)-1
            bnum = Tdomain%specel(pelems(n))%lnum/VCHUNK
            ee = mod(Tdomain%specel(pelems(n))%lnum, VCHUNK)
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        dx = c * (1d0/Tdomain%spmldom%PMLDumpSx_(i,j,k,1,bnum,ee) - 1d0)
                        dy = c * (1d0/Tdomain%spmldom%PMLDumpSy_(i,j,k,1,bnum,ee) - 1d0)
                        dz = c * (1d0/Tdomain%spmldom%PMLDumpSz_(i,j,k,1,bnum,ee) - 1d0)
                        acoef(0,ip) = dx + dy + dz
                        acoef(1,ip) = dx*dy + dx*dz + dy*dz
                        acoef(2,ip) = dx*dy*dz
                        ip = ip + 1
                    end do
                end do
            end do
        end do
#endif
    end subroutine pml_coefs

    !> Deformations aux points des nel elements solides ecrits (elb, ele) du champ U, dans
    !! out(:,0:nel*ngll**3-1), meme ordre et memes operations que physical_part_deriv_ijk.
    subroutine strain_sol(dom, nel, U, out)
        type(domain_solid), intent(in) :: dom
        integer, intent(in) :: nel
        real(fpp), dimension(0:,0:), intent(in) :: U
        real(kind=8), dimension(0:,0:), intent(inout) :: out
        !
        integer :: n, i, j, k, l, c, ngll, bnum, ee, ip
        real(fpp) :: dxi, deta, dzeta, g(0:2,0:2)

        ngll = dom%ngll
        !$acc parallel loop collapse(4) async(1) firstprivate(ngll) &
        !$acc& private(bnum, ee, ip, l, c, dxi, deta, dzeta, g) &
        !$acc& present(dom, dom%m_Idom, dom%m_InvGrad, dom%hprime, elb, ele, U, out)
        do n = 0, nel-1
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        bnum = elb(n)
                        ee = ele(n)
                        ip = ((n*ngll + k)*ngll + j)*ngll + i
                        ! g(a,c) = dU_c/dA
                        do c = 0, 2
                            dxi = 0d0
                            deta = 0d0
                            dzeta = 0d0
                            do l = 0, ngll-1
                                dxi = dxi + U(dom%Idom_(l,j,k,bnum,ee),c)*dom%hprime(l,i)
                                deta = deta + U(dom%Idom_(i,l,k,bnum,ee),c)*dom%hprime(l,j)
                                dzeta = dzeta + U(dom%Idom_(i,j,l,bnum,ee),c)*dom%hprime(l,k)
                            end do
                            g(0,c) = dxi*dom%InvGrad_(0,0,i,j,k,bnum,ee) + deta*dom%InvGrad_(0,1,i,j,k,bnum,ee) &
                                + dzeta*dom%InvGrad_(0,2,i,j,k,bnum,ee)
                            g(1,c) = dxi*dom%InvGrad_(1,0,i,j,k,bnum,ee) + deta*dom%InvGrad_(1,1,i,j,k,bnum,ee) &
                                + dzeta*dom%InvGrad_(1,2,i,j,k,bnum,ee)
                            g(2,c) = dxi*dom%InvGrad_(2,0,i,j,k,bnum,ee) + deta*dom%InvGrad_(2,1,i,j,k,bnum,ee) &
                                + dzeta*dom%InvGrad_(2,2,i,j,k,bnum,ee)
                        end do
                        out(0,ip) = g(0,0)
                        out(1,ip) = g(1,1)
                        out(2,ip) = g(2,2)
                        out(3,ip) = 0.5d0*(g(1,0)+g(0,1))
                        out(4,ip) = 0.5d0*(g(2,0)+g(0,2))
                        out(5,ip) = 0.5d0*(g(2,1)+g(1,2))
                    end do
                end do
            end do
        end do
    end subroutine strain_sol

    !> Appele a chaque pas de temps, apres le calcul des champs (file OpenACC 1).
    subroutine save_strain_output(Tdomain)
        type(domain), intent(inout) :: Tdomain
        !
        integer :: ip, r, kf, c, npp, nb
        integer(kind=8) :: pos
        real(kind=8) :: dtp, tu, tw

        if (.not. strain_open) return
        nstep = nstep + 1
        if (mod(nstep, period) /= 0) return
        npp = npts - npts_s
        dtp = Tdomain%TimeD%dtmin * period
        select case (mode)
        case (0)
            call strain_of(Tdomain, Tdomain%sdom%champs(0)%Depla, buffer)
            nb = npts
            !$acc update host(buffer(:,0:nb-1)) async(1) if(nb > 0)
            !$acc wait(1)
            if (npts > 0) write(strain_unit) buffer(:,0:npts-1)
        case (1)
            call gather_nodes(Tdomain%sdom%champs(0)%Depla)
            if (npp > 0) then
                call strain_pml_of(Tdomain, buffer)
                !$acc parallel loop async(1) firstprivate(dtp) private(c) present(buffer, iint, acoef, epf)
                do ip = 0, npp-1
                    do c = 0, 5
                        iint(c,ip,1) = iint(c,ip,1) + dtp * buffer(c,npts_s+ip)
                        iint(c,ip,2) = iint(c,ip,2) + dtp * iint(c,ip,1)
                        iint(c,ip,3) = iint(c,ip,3) + dtp * iint(c,ip,2)
                        epf(c,ip) = real(buffer(c,npts_s+ip) + acoef(0,ip)*iint(c,ip,1) &
                            + acoef(1,ip)*iint(c,ip,2) + acoef(2,ip)*iint(c,ip,3), 4)
                    end do
                end do
                !$acc update host(epf) async(1)
            end if
            !$acc wait(1)
            if (nnod > 0) write(disp_unit) dispf(:,0:nnod-1)
            if (npp > 0) write(epml_unit) epf
        case (2)
            r = nstep/period - 1
            kf = nrecf - 2 - r
            if (kf < 0 .or. npts == 0) return
            call strain_of(Tdomain, Tdomain%sdom%champs(0)%Depla, buffer)
            if (nnod > 0) then
                pos = int(kf, 8) * int(nnod, 8) * 24_8 + 1_8
                read(disp_unit, pos=pos) dispf(:,0:nnod-1)
                call scatter_nodes()
            end if
            call strain_of(Tdomain, ufull, bufferf)
            if (npp > 0) then
                call strain_pml_of(Tdomain, buffer)
                pos = int(kf, 8) * int(npp, 8) * 24_8 + 1_8
                read(epml_unit, pos=pos) epf
                !$acc update device(epf) async(1)
                !$acc parallel loop async(1) private(c) present(bufferf, epf)
                do ip = 0, npp-1
                    do c = 0, 5
                        bufferf(c,npts_s+ip) = real(epf(c,ip), 8)
                    end do
                end do
            end if
            !$acc parallel loop async(1) firstprivate(dtp) private(tu, tw) present(buffer, bufferf, kl, km)
            do ip = 0, npts-1
                tu = bufferf(0,ip) + bufferf(1,ip) + bufferf(2,ip)
                tw = buffer(0,ip) + buffer(1,ip) + buffer(2,ip)
                kl(ip) = kl(ip) - dtp * tu * tw
                km(ip) = km(ip) - dtp * 2d0 * (bufferf(0,ip)*buffer(0,ip) + bufferf(1,ip)*buffer(1,ip) &
                    + bufferf(2,ip)*buffer(2,ip) + 2d0*(bufferf(3,ip)*buffer(3,ip) &
                    + bufferf(4,ip)*buffer(4,ip) + bufferf(5,ip)*buffer(5,ip)))
            end do
            ! dispf et epf sont relus sur l''hote au prochain enregistrement
            !$acc wait(1)
        end select
    end subroutine save_strain_output

    !> dispf = U aux noeuds ecrits, rapatrie sur l''hote (file 1, a attendre).
    subroutine gather_nodes(U)
        real(fpp), dimension(0:,0:), intent(in) :: U
        integer :: m, n
        n = nnod
        if (n == 0) return
        !$acc parallel loop async(1) present(nodes, U, dispf)
        do m = 0, n-1
            dispf(0,m) = U(nodes(m),0)
            dispf(1,m) = U(nodes(m),1)
            dispf(2,m) = U(nodes(m),2)
        end do
        !$acc update host(dispf(:,0:n-1)) async(1)
    end subroutine gather_nodes

    !> ufull = dispf aux noeuds ecrits (dispf lu sur l''hote).
    subroutine scatter_nodes()
        integer :: m, n
        n = nnod
        !$acc update device(dispf(:,0:n-1)) async(1)
        !$acc parallel loop async(1) present(nodes, ufull, dispf)
        do m = 0, n-1
            ufull(nodes(m),0) = dispf(0,m)
            ufull(nodes(m),1) = dispf(1,m)
            ufull(nodes(m),2) = dispf(2,m)
        end do
    end subroutine scatter_nodes

    subroutine close_strain_output()
        integer :: ku
        if (strain_open) then
            if (mode == 0) then
                close(strain_unit)
            else
                close(disp_unit)
                if (npts > npts_s) close(epml_unit)
            end if
            if (mode == 2 .and. npts > 0) then
                !$acc update host(kl, km)
                open(newunit=ku, file="strain/kernel."//trim(adjustl(strrank(rank_of_kernels)))//".bin", &
                    access="stream", form="unformatted", status="replace", action="write")
                write(ku) kl(0:npts-1), km(0:npts-1)
                close(ku)
            end if
            if (mode /= 0) then
                if (mode == 2) then
                    !$acc exit data delete(kl, km, ufull, bufferf)
                end if
                if (npts > npts_s) then
                    if (mode == 1) then
                        !$acc exit data delete(acoef, iint)
                    end if
                    !$acc exit data delete(epf)
                end if
                !$acc exit data delete(nodes, dispf)
            end if
            !$acc exit data delete(elb, ele, pelb, pele, buffer)
            strain_open = .false.
        end if
        if (allocated(elb)) deallocate(elb, ele, pelb, pele)
        if (allocated(nodes)) deallocate(nodes)
        if (allocated(nodemap)) deallocate(nodemap)
        if (allocated(dispf)) deallocate(dispf)
        if (allocated(bufferf)) deallocate(bufferf)
        if (allocated(ufull)) deallocate(ufull)
        if (allocated(kl)) deallocate(kl)
        if (allocated(km)) deallocate(km)
        if (allocated(buffer)) deallocate(buffer)
        if (allocated(elems)) deallocate(elems)
        if (allocated(pelems)) deallocate(pelems)
        if (allocated(epf)) deallocate(epf)
        if (allocated(acoef)) deallocate(acoef)
        if (allocated(iint)) deallocate(iint)
    end subroutine close_strain_output

end module mstrain_output

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
