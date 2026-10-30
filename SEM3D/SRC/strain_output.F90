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
    integer, allocatable, dimension(:), save :: elems
    real(kind=8), allocatable, dimension(:,:), save :: buffer
    integer, save :: mode = 0, nnod = 0, nrecf = 0, disp_unit, rank_of_kernels = 0
    integer, allocatable, dimension(:), save :: nodes      ! noeuds du domaine solide ecrits
    integer, allocatable, dimension(:), save :: nodemap    ! indice du domaine -> rang dans nodes
    real(kind=8), allocatable, dimension(:,:), save :: dispf, bufferf, ufull
    real(kind=8), allocatable, dimension(:), save :: kl, km

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
        character(len=MAX_FILE_SIZE) :: fgeom, fdata
        real(fpp), allocatable, dimension(:,:) :: coords
        real(fpp), allocatable, dimension(:) :: weight
        integer, allocatable, dimension(:) :: elem_of_pt
        real(fpp), dimension(:), allocatable :: tmp

        if (Tdomain%rank == 0) then
            ierr = sem_mkdir("strain")
        end if
        call MPI_Barrier(Tdomain%communicateur, ierr)

        nel = 0
        do n = 0, Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain == DM_SOLID_CG) nel = nel + 1
        end do
        allocate(elems(0:nel-1))
        nel = 0
        do n = 0, Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain == DM_SOLID_CG) then
                elems(nel) = n
                nel = nel + 1
            end if
        end do

        ngll = Tdomain%sdom%ngll
        npts = nel*ngll*ngll*ngll
        allocate(buffer(0:5, 0:max(npts,1)-1))
        allocate(coords(0:2, 0:max(npts,1)-1), weight(0:max(npts,1)-1), elem_of_pt(0:max(npts,1)-1))

        ip = 0
        do n = 0, nel-1
            bnum = Tdomain%specel(elems(n))%lnum/VCHUNK
            ee = mod(Tdomain%specel(elems(n))%lnum, VCHUNK)
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        coords(:,ip) = Tdomain%GlobCoord(:, Tdomain%specel(elems(n))%Iglobnum(i,j,k))
                        weight(ip) = Tdomain%sdom%GLLw(i)*Tdomain%sdom%GLLw(j)*Tdomain%sdom%GLLw(k) &
                            * Tdomain%sdom%Jacob_(i,j,k,bnum,ee)
                        elem_of_pt(ip) = elems(n)
                        ip = ip + 1
                    end do
                end do
            end do
        end do

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
        period = max(1, Tdomain%config%strain_period)
        nstep = 0
        call write_attr_int(fid, "period", period)
        call h5fclose_f(fid, hdferr)
        deallocate(coords, weight, elem_of_pt)

        mode = Tdomain%config%strain_mode
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
        end if
        strain_open = .true.
        if (Tdomain%rank == 0) write(*,*) "--> STRAIN OUTPUT: ", nel, " elements, ", npts, " points on rank 0"
    end subroutine init_strain_output

    subroutine init_nodes(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: n, i, j, k, ngll, bnum, ee, ind
        ngll = Tdomain%sdom%ngll
        allocate(nodemap(0:Tdomain%sdom%nglltot))
        nodemap = -1
        nnod = 0
        do n = 0, size(elems)-1
            bnum = Tdomain%specel(elems(n))%lnum/VCHUNK
            ee = mod(Tdomain%specel(elems(n))%lnum, VCHUNK)
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        ind = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                        if (nodemap(ind) < 0) then
                            nodemap(ind) = nnod
                            nnod = nnod + 1
                        end if
                    end do
                end do
            end do
        end do
        allocate(nodes(0:max(nnod,1)-1))
        do ind = 0, Tdomain%sdom%nglltot
            if (nodemap(ind) >= 0) nodes(nodemap(ind)) = ind
        end do
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

    !> Deformations aux points ecrits du champ de deplacement U (indice du domaine solide).
    subroutine strain_of(Tdomain, U, out)
        use deriv3d
        type(domain), intent(inout) :: Tdomain
        real(fpp), dimension(0:,0:), intent(in) :: U
        real(kind=8), dimension(0:,0:), intent(out) :: out
        !
        integer :: n, i, j, k, ip, ngll, bnum, ee, ind
        real(fpp), dimension(0:Tdomain%sdom%ngll-1,0:Tdomain%sdom%ngll-1,0:Tdomain%sdom%ngll-1,0:2) :: fieldU
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        real(fpp) :: DXX, DXY, DXZ, DYX, DYY, DYZ, DZX, DZY, DZZ

        ngll = Tdomain%sdom%ngll
        ip = 0
        do n = 0, size(elems)-1
            bnum = Tdomain%specel(elems(n))%lnum/VCHUNK
            ee = mod(Tdomain%specel(elems(n))%lnum, VCHUNK)
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        ind = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                        fieldU(i,j,k,:) = U(ind,:)
                    end do
                end do
            end do
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        ! DAB = dU_B/dA
                        invgrad_ijk = Tdomain%sdom%InvGrad_(:,:,i,j,k,bnum,ee)
                        call physical_part_deriv_ijk(i,j,k,ngll,Tdomain%sdom%hprime,&
                            invgrad_ijk,fieldU(:,:,:,0),DXX,DYX,DZX)
                        call physical_part_deriv_ijk(i,j,k,ngll,Tdomain%sdom%hprime,&
                            invgrad_ijk,fieldU(:,:,:,1),DXY,DYY,DZY)
                        call physical_part_deriv_ijk(i,j,k,ngll,Tdomain%sdom%hprime,&
                            invgrad_ijk,fieldU(:,:,:,2),DXZ,DYZ,DZZ)
                        out(0,ip) = DXX
                        out(1,ip) = DYY
                        out(2,ip) = DZZ
                        out(3,ip) = 0.5d0*(DYX+DXY)
                        out(4,ip) = 0.5d0*(DZX+DXZ)
                        out(5,ip) = 0.5d0*(DZY+DYZ)
                        ip = ip + 1
                    end do
                end do
            end do
        end do
    end subroutine strain_of

    subroutine save_strain_output(Tdomain)
        type(domain), intent(inout) :: Tdomain
        !
        integer :: m, ip, r, kf
        integer(kind=8) :: pos
        real(kind=8) :: dtp, tu, tw

        if (.not. strain_open) return
        nstep = nstep + 1
        if (mod(nstep, period) /= 0) return
        !$acc update host(Tdomain%sdom%champs(0)%Depla) async(2) wait(1)
        !$acc wait(2)
        select case (mode)
        case (0)
            call strain_of(Tdomain, Tdomain%sdom%champs(0)%Depla, buffer)
            if (npts > 0) write(strain_unit) buffer(:,0:npts-1)
        case (1)
            do m = 0, nnod-1
                dispf(:,m) = Tdomain%sdom%champs(0)%Depla(nodes(m),:)
            end do
            if (nnod > 0) write(disp_unit) dispf(:,0:nnod-1)
        case (2)
            r = nstep/period - 1
            kf = nrecf - 2 - r
            if (kf < 0 .or. npts == 0) return
            call strain_of(Tdomain, Tdomain%sdom%champs(0)%Depla, buffer)
            pos = int(kf, 8) * int(nnod, 8) * 24_8 + 1_8
            read(disp_unit, pos=pos) dispf(:,0:nnod-1)
            do m = 0, nnod-1
                ufull(nodes(m),:) = dispf(:,m)
            end do
            call strain_of(Tdomain, ufull, bufferf)
            dtp = Tdomain%TimeD%dtmin * period
            do ip = 0, npts-1
                tu = bufferf(0,ip) + bufferf(1,ip) + bufferf(2,ip)
                tw = buffer(0,ip) + buffer(1,ip) + buffer(2,ip)
                kl(ip) = kl(ip) - dtp * tu * tw
                km(ip) = km(ip) - dtp * 2d0 * (bufferf(0,ip)*buffer(0,ip) + bufferf(1,ip)*buffer(1,ip) &
                    + bufferf(2,ip)*buffer(2,ip) + 2d0*(bufferf(3,ip)*buffer(3,ip) &
                    + bufferf(4,ip)*buffer(4,ip) + bufferf(5,ip)*buffer(5,ip)))
            end do
        end select
    end subroutine save_strain_output

    subroutine close_strain_output()
        integer :: ku
        if (strain_open) then
            if (mode == 0) then
                close(strain_unit)
            else
                close(disp_unit)
            end if
            if (mode == 2 .and. npts > 0) then
                open(newunit=ku, file="strain/kernel."//trim(adjustl(strrank(rank_of_kernels)))//".bin", &
                    access="stream", form="unformatted", status="replace", action="write")
                write(ku) kl(0:npts-1), km(0:npts-1)
                close(ku)
            end if
            strain_open = .false.
        end if
        if (allocated(nodes)) deallocate(nodes)
        if (allocated(nodemap)) deallocate(nodemap)
        if (allocated(dispf)) deallocate(dispf)
        if (allocated(bufferf)) deallocate(bufferf)
        if (allocated(ufull)) deallocate(ufull)
        if (allocated(kl)) deallocate(kl)
        if (allocated(km)) deallocate(km)
        if (allocated(buffer)) deallocate(buffer)
        if (allocated(elems)) deallocate(elems)
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
