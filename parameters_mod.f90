! filepath: /home/yanel/PA/fortran/first/exercices/heat/parameters_mod.f90
! -----------------------------------------------------------------------
! Module de paramètres pour la résolution de l'équation de la chaleur
! -----------------------------------------------------------------------

module parameters_mod
  implicit none
  
  ! Précision
  integer, parameter :: dp = kind(1.0d0)
  
  ! Paramètres du problème
  integer, parameter :: N = 100                      ! Nombre de points de grille (total)
  real(dp), parameter :: L = 1.0_dp                  ! Longueur du domaine
  real(dp), parameter :: alpha = 1.0_dp              ! Coefficient de diffusion
  real(dp), parameter :: dt = 0.001_dp               ! Pas de temps
  real(dp), parameter :: t_final = 1.0_dp            ! Temps final
  real(dp), parameter :: epsilon = 1.0e-6_dp         ! Tolérance pour la convergence CG
  integer, parameter :: max_iter_cg = 1000           ! Nombre maximum d'itérations pour CG
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)  ! π
  
  ! Paramètres dérivés (calculés lors de la compilation)
  real(dp), parameter :: dx = L / real(N, dp)
  real(dp), parameter :: lambda = alpha * dt / (dx * dx)
  
end module parameters_mod