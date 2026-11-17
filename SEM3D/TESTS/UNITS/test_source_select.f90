!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!! Tests du choix de l'element qui porte une source (msource_select).
!! Le programme s'arrete avec le code 1 au premier test en echec.
program test_source_select

    use constants
    use msource_select
    implicit none

    real(fpp), parameter :: EPS = 1d-10
    real(fpp), dimension(0:2,4) :: lc
    real(fpp), dimension(0:2,4) :: ce
    integer, dimension(4) :: doms
    integer :: p(4), q

    ! Centres des elements candidats, distincts (ils ne departagent que des candidats valables)
    ce(:,1) = (/ 0d0, 0d0, 0d0 /)
    ce(:,2) = (/ 1d0, 0d0, 0d0 /)
    ce(:,3) = (/ 0d0, 1d0, 0d0 /)
    ce(:,4) = (/ 1d0, 1d0, 0d0 /)

    ! Un seul candidat, a l'interieur, du bon domaine
    lc(:,1) = (/ 0.2d0, -0.3d0, 0.5d0 /)
    doms(1) = DM_SOLID_CG
    call check("interieur, source solide", select_source_candidate(1, lc, ce, doms, .false., EPS), 1)
    call check("interieur, source fluide dans un solide", select_source_candidate(1, lc, ce, doms, .true., EPS), 0)

    ! Point sur une interface fluide-solide : l'element fluide est trouve en premier
    lc(:,1) = (/ 0.2d0, -0.3d0, -1d0 /)
    lc(:,2) = (/ 0.2d0, -0.3d0,  1d0 /)
    doms(1:2) = (/ DM_FLUID_CG, DM_SOLID_CG /)
    call check("interface, source solide", select_source_candidate(2, lc, ce, doms, .false., EPS), 2)
    call check("interface, source fluide", select_source_candidate(2, lc, ce, doms, .true., EPS), 1)
    ! Ordre inverse
    doms(1:2) = (/ DM_SOLID_CG, DM_FLUID_CG /)
    call check("interface inversee, source solide", select_source_candidate(2, lc, ce, doms, .false., EPS), 1)
    call check("interface inversee, source fluide", select_source_candidate(2, lc, ce, doms, .true., EPS), 2)

    ! Domaines PML et DG
    doms(1:2) = (/ DM_FLUID_CG_PML, DM_SOLID_CG_PML /)
    call check("interface PML, source solide", select_source_candidate(2, lc, ce, doms, .false., EPS), 2)
    doms(1:2) = (/ DM_SOLID_DG, DM_FLUID_DG /)
    call check("interface DG, source fluide", select_source_candidate(2, lc, ce, doms, .true., EPS), 2)

    ! Candidat hors de l'element (coordonnees locales hors de [-1,1]) puis candidat valable
    lc(:,1) = (/ 1.5d0, 0d0, 0d0 /)
    lc(:,2) = (/ 0.1d0, 0d0, 0d0 /)
    doms(1:2) = DM_SOLID_CG
    call check("premier candidat hors de l'element", select_source_candidate(2, lc, ce, doms, .false., EPS), 2)

    ! Tolerance : sommet a 1e-12 pres
    lc(:,1) = (/ 1d0 + 1d-12, -1d0 - 1d-12, 1d0 /)
    call check("sommet, dans la tolerance", select_source_candidate(1, lc, ce, doms, .false., EPS), 1)
    lc(:,1) = (/ 1d0 + 1d-8, 0d0, 0d0 /)
    call check("au-dela de la tolerance", select_source_candidate(1, lc, ce, doms, .false., EPS), 0)

    ! Sommet partage par 4 elements, seul le dernier est solide
    lc(:,1) = (/ -1d0, -1d0, 1d0 /)
    lc(:,2) = (/  1d0, -1d0, 1d0 /)
    lc(:,3) = (/ -1d0,  1d0, 1d0 /)
    lc(:,4) = (/  1d0,  1d0, 1d0 /)
    doms = (/ DM_FLUID_CG, DM_FLUID_CG, DM_FLUID_CG_PML, DM_SOLID_CG /)
    call check("sommet, dernier candidat solide", select_source_candidate(4, lc, ce, doms, .false., EPS), 4)
    doms(4) = DM_FLUID_CG
    call check("sommet, aucun candidat solide", select_source_candidate(4, lc, ce, doms, .false., EPS), 0)

    ! Aucun candidat
    call check("aucun candidat", select_source_candidate(0, lc, ce, doms, .false., EPS), 0)

    ! Plusieurs candidats valables : le choix ne depend pas de leur ordre (decoupage MPI)
    ! Face horizontale : l'element du dessus
    lc(:,1) = (/ 0.2d0, -0.3d0, -1d0 /)
    lc(:,2) = (/ 0.2d0, -0.3d0,  1d0 /)
    ce(:,1) = (/ 0d0, 0d0,  50d0 /)
    ce(:,2) = (/ 0d0, 0d0, -50d0 /)
    doms(1:2) = DM_SOLID_CG
    call check("face horizontale, dessus en premier", select_source_candidate(2, lc, ce, doms, .false., EPS), 1)
    call check("face horizontale, dessus en second", &
        select_source_candidate(2, lc(:,(/2,1/)), ce(:,(/2,1/)), doms, .false., EPS), 2)
    ! La source la plus interieure l'emporte sur l'element le plus haut
    lc(:,1) = (/ 1d0, 0d0, 0d0 /)
    lc(:,2) = (/ 0.5d0, 0d0, 0d0 /)
    call check("plus interieure", select_source_candidate(2, lc, ce, doms, .false., EPS), 2)
    ! Arete verticale partagee par 4 elements a la meme hauteur : le plus petit x, puis y
    lc(:,1) = (/  1d0,  1d0, 0d0 /)
    lc(:,2) = (/ -1d0,  1d0, 0d0 /)
    lc(:,3) = (/  1d0, -1d0, 0d0 /)
    lc(:,4) = (/ -1d0, -1d0, 0d0 /)
    ce(:,1) = (/ -50d0, -50d0, 0d0 /)
    ce(:,2) = (/  50d0, -50d0, 0d0 /)
    ce(:,3) = (/ -50d0,  50d0, 0d0 /)
    ce(:,4) = (/  50d0,  50d0, 0d0 /)
    doms = DM_SOLID_CG
    ! Toutes les permutations de 4 candidats donnent l'element 1 (x puis y les plus petits)
    do q = 0, 23
        p = permutation(q)
        call check("arete, permutation", p(select_source_candidate(4, lc(:,p), ce(:,p), doms(p), .false., EPS)), 1)
    end do

    write(*,*) "test_source_select: OK"

contains

    !> q-ieme permutation de (1, 2, 3, 4), q dans 0..23
    function permutation(q) result(p)
        integer, intent(in) :: q
        integer :: p(4), rest(4), r, k, n
        rest = (/ 1, 2, 3, 4 /)
        r = q
        n = 4
        do k = 1, 4
            p(k) = rest(mod(r, n) + 1)
            rest(mod(r, n) + 1:n - 1) = rest(mod(r, n) + 2:n)
            r = r / n
            n = n - 1
        end do
    end function permutation

    subroutine check(name, got, expected)
        character(len=*), intent(in) :: name
        integer, intent(in) :: got, expected
        if (got /= expected) then
            write(*,*) "ECHEC ", name, " : obtenu ", got, ", attendu ", expected
            stop 1
        end if
    end subroutine check

end program test_source_select
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
