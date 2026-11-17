!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file source_select.f90
!! \brief Choix de l''element qui porte une source ponctuelle
!!
!<

module msource_select
    use constants
    implicit none
contains

    !> Vrai pour un domaine fluide (meme classement que les sources : tout ce qui n''est pas solide)
    pure logical function domain_is_fluid(dom)
        integer, intent(in) :: dom
        domain_is_fluid = .not. (dom == DM_SOLID_CG .or. dom == DM_SOLID_DG .or. dom == DM_SOLID_CG_PML)
    end function domain_is_fluid

    !> Cle de choix d'un candidat : (max |xi|, |eta|, |zeta|, -z, x, y du centre de l'element).
    !! La plus petite cle l'emporte : l'element ou la source est la plus interieure, puis le plus
    !! haut, puis le plus petit en x, puis en y. Elle ne depend que de la geometrie, donc pas
    !! du decoupage MPI ni de la numerotation locale des elements.
    pure function source_candidate_key(localcoord, centre) result(key)
        real(fpp), dimension(0:2), intent(in) :: localcoord, centre
        real(8), dimension(0:3) :: key
        key(0) = maxval(abs(localcoord))
        key(1) = -centre(2)
        key(2) = centre(0)
        key(3) = centre(1)
    end function source_candidate_key

    !> a < b pour les cles de source_candidate_key : la premiere composante a 1e-8 pres
    !! (coordonnees locales), les suivantes a eps pres (coordonnees physiques)
    pure logical function source_key_less(a, b, eps)
        real(8), dimension(0:3), intent(in) :: a, b
        real(fpp), intent(in) :: eps
        real(8) :: tol
        integer :: k
        source_key_less = .false.
        do k = 0, 3
            tol = eps
            if (k == 0) tol = 1d-8
            if (a(k) < b(k) - tol) then
                source_key_less = .true.
                return
            else if (a(k) > b(k) + tol) then
                return
            end if
        end do
    end function source_key_less

    !> Indice (1..n) du candidat qui contient le point (coordonnees locales dans [-1-eps, 1+eps]),
    !! dont le domaine convient au type de source, et de plus petite cle (source_candidate_key,
    !! centre(:,i) le centre de l'element i) ; 0 si aucun.
    !! Un point sur une frontiere a plusieurs candidats : on ecarte ceux dont le domaine ne
    !! convient pas (source solide sur une interface fluide-solide par exemple), puis la cle
    !! choisit le meme element quel que soit l'ordre des candidats. Les gradients des fonctions
    !! de forme d'une source en tenseur de moment different d'un element a l'autre, et une source
    !! sur une interface prend le materiau de l'element choisi.
    pure integer function select_source_candidate(n, localcoord, centre, doms, fluid_source, eps)
        integer, intent(in) :: n
        real(fpp), dimension(0:2,n), intent(in) :: localcoord, centre
        integer, dimension(n), intent(in) :: doms
        logical, intent(in) :: fluid_source
        real(fpp), intent(in) :: eps
        real(8), dimension(0:3) :: key, best
        integer :: i

        select_source_candidate = 0
        do i = 1,n
            if (any(localcoord(:,i) < -1-eps) .or. any(localcoord(:,i) > 1+eps)) cycle
            if (domain_is_fluid(doms(i)) .neqv. fluid_source) cycle
            key = source_candidate_key(localcoord(:,i), centre(:,i))
            if (select_source_candidate == 0) then
                best = key
                select_source_candidate = i
            else if (source_key_less(key, best, eps)) then
                best = key
                select_source_candidate = i
            end if
        end do
    end function select_source_candidate

end module msource_select
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
