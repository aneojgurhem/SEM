!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file SourcePosition.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine SourcePosition (Tdomain)
    use sdomain
    use ssources
    use sem_mpi
    use mshape8
    use mshape27
    use mlinalg
    use mlocations3d
    use msource_select
    implicit none
    type (domain), intent(inout) :: Tdomain
    !
    integer :: rg, ierr, nnodes
    integer :: n_src, i,j,k
    real(fpp) :: xs,ys,zs,xi,eta,zeta
    real(fpp), dimension(0:2,0:2) :: LocInvGrad
    real(fpp), dimension(:,:), allocatable :: coord
    real(fpp) :: R
    integer :: src_proc, nmax, n_el, nb_src_inproc, attempt
    integer, parameter :: NMAXEL=20
    integer, dimension(NMAXEL) :: elems, doms
    real(fpp), dimension(0:2,NMAXEL) :: coordloc, centre
    real(8), dimension(0:3) :: key
    real(8), dimension(:,:), allocatable :: keys
    real(fpp), parameter :: EPSD = 1D-10
    real(fpp) :: EPS
    type(source), dimension(:), allocatable :: ssource_temp
    logical, dimension(:), allocatable :: src_inproc

    rg = Tdomain%rank
    nb_src_inproc = 0
    nnodes = Tdomain%n_nodes
    allocate(coord(0:2, 0:nnodes-1))
    allocate(src_inproc(0:Tdomain%n_source-1))
    src_inproc = .false.
    EPS=EPSD*Tdomain%dxmax
    do n_src = 0, Tdomain%n_source-1
        xs = Tdomain%Ssource(n_src)%Xsource
        ys = Tdomain%Ssource(n_src)%Ysource
        zs = Tdomain%Ssource(n_src)%Zsource

        ! On trouve les elements autour de la source. Si aucun processeur ne la trouve (source sur
        ! un noeud de la surface, au-dessus des boites englobantes d'une erreur d'arrondi), on
        ! cherche parmi les elements autour du noeud le plus proche.
        do attempt = 1, 2
        nmax = NMAXEL
        if (attempt == 1) then
            call find_location(Tdomain, xs, ys, zs, nmax, elems, coordloc, n_src)
        else
            call find_location_closest(Tdomain, xs, ys, zs, nmax, elems, coordloc, n_src)
        end if
        n_el = -1

        ! Element qui contient la source, dont le domaine convient (source fluide dans un
        ! element fluide, sinon solide) et de plus petite cle (source_candidate_key). Un point
        ! sur une frontiere a plusieurs candidats : on ne s''arrete pas au premier.
        do i=1,nmax
            doms(i) = Tdomain%specel(elems(i))%domain
            do j = 0, nnodes-1
                coord(0:2, j) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(elems(i))%Control_Nodes(j))
            enddo
            centre(:,i) = sum(coord, dim=2) / nnodes
        end do
        k = select_source_candidate(nmax, coordloc(:,1:nmax), centre(:,1:nmax), doms(1:nmax), &
            Tdomain%sSource(n_src)%i_type_source == 3, EPS)
        key = huge(1d0)
        if (k > 0) then
            n_el = elems(k)
            xi   = coordloc(0,k)
            eta  = coordloc(1,k)
            zeta = coordloc(2,k)
            key = source_candidate_key(coordloc(:,k), centre(:,k))
        end if

        Tdomain%Ssource(n_src)%elem = n_el
        ! On trouve le processeur qui contient la src : celui de plus petite cle, la meme
        ! sur tous les processeurs, donc le meme element quel que soit le decoupage
        allocate(keys(0:3, 0:Tdomain%nb_procs-1))
        call MPI_Allgather(key, 4, MPI_DOUBLE_PRECISION, keys, 4, MPI_DOUBLE_PRECISION, &
                           Tdomain%communicateur, ierr)
        src_proc = -1
        do j = 0, Tdomain%nb_procs-1
            if (keys(0,j) == huge(1d0)) cycle
            if (src_proc == -1) then
                src_proc = j
            else if (source_key_less(keys(:,j), keys(:,src_proc), EPS)) then
                src_proc = j
            end if
        end do
        deallocate(keys)
        if (src_proc /= -1) exit
        end do
        if (rg==src_proc) then
            Tdomain%sSource(n_src)%proc = rg
        else
            Tdomain%sSource(n_src)%proc = -1
        end if

        if (src_proc==-1) then
            if (rg==0) then
                write(*,*) "The source at (", xs, ys, zs, ") doesn't appear to be on any processor"
                write(*,*) "Please verify that the source location is within the computation domain"
                write(*,*) "And the source type (solid, or fluid) is in accordance with its location"

                write(*,*) "ALSO PRINTING n_el and proc   ", n_el, Tdomain%sSource(n_src)%proc

            end if
            stop 1
        endif

        if (rg==src_proc) then
            n_el = Tdomain%sSource(n_src)%elem
            do j = 0, nnodes-1
                coord(0:2, j) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(n_el)%Control_Nodes(j))
            enddo
            if (nnodes==27) then
                call shape27_local2jacob(coord, xi, eta, zeta, LocInvGrad)
            else if (nnodes==8) then
                call shape8_local2jacob(coord, xi, eta, zeta, LocInvGrad)
            endif
            Tdomain%sSource(n_src)%refcoord(0) = xi
            Tdomain%sSource(n_src)%refcoord(1) = eta
            Tdomain%sSource(n_src)%refcoord(2) = zeta
            if (Tdomain%config%verbose_level .gt. 1) then
                write(*,*) "Found source", n_src, "on proc", src_proc
                write(*,*) "Source local coordinates:", Tdomain%sSource(n_src)%refcoord
                write(*,*) "From position", xs, ys, zs
            end if
            call invert_3d (LocInvGrad, R)
            Tdomain%sSource(n_src)%InvGrad(0:2,0:2) = LocInvGrad(0:2,0:2)
            src_inproc(n_src) = .true.
            nb_src_inproc = nb_src_inproc +1
        endif
    enddo

    ! Reallocation des vecteurs de sources ne comprenant que les sources
    ! reellement presentes dans le domaine traitee par le processeur
    i = 0
    if (nb_src_inproc>0) then
        allocate(ssource_temp(0:nb_src_inproc-1))
        do n_src = 0, Tdomain%n_source-1
            if (src_inproc(n_src)) then
                ssource_temp(i) = Tdomain%sSource(n_src)
                i = i+1
            endif
        enddo
        deallocate(Tdomain%sSource)
        allocate(Tdomain%sSource(0:nb_src_inproc-1))
        Tdomain%sSource = ssource_temp
        deallocate(ssource_temp)
    else ! No source in current proc
        Tdomain%logicD%any_source = .false.
        deallocate(Tdomain%sSource)
    endif
    Tdomain%n_source = nb_src_inproc

end subroutine SourcePosition

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
