! -----------------------------------------------------------------------
! Module de méthodes numériques pour résoudre l'équation de la chaleur
! -----------------------------------------------------------------------

module solver_mod
  use mpi
  use parameters_mod, only: dp, lambda, max_iter_cg, epsilon
  implicit none
  private
  
  ! Procédures exportées
  public :: implicit_step, conjugate_gradient
  
contains

  ! Effectue un pas de temps implicite
  subroutine implicit_step(u_loc, u_new_loc, N_loc, alpha, dt, myrank, nprocs, &
                          left_neighbor, right_neighbor)
    integer, intent(in) :: N_loc, myrank, nprocs, left_neighbor, right_neighbor
    real(dp), intent(in) :: u_loc(0:N_loc+1), alpha, dt
    real(dp), intent(out) :: u_new_loc(0:N_loc+1)
    real(dp) :: a, b, c
    real(dp), allocatable :: b_loc(:), x_loc(:)
    
    ! Coefficients de la matrice tridiagonale
    a = -lambda  ! Coefficient sous-diagonal
    b = 1.0_dp + 2.0_dp * lambda  ! Coefficient diagonal
    c = -lambda  ! Coefficient sur-diagonal
    
    ! Allouer les vecteurs pour le système linéaire
    allocate(b_loc(N_loc), x_loc(N_loc))
    
    ! Construction du vecteur b local (second membre)
    b_loc(:) = u_loc(1:N_loc)
    
    ! Appliquer les conditions aux limites au second membre
    if (myrank == 0) then
      b_loc(1) = b_loc(1) + lambda * 0.0_dp  ! u(0,t) = 0
    end if
    if (myrank == nprocs-1) then
      b_loc(N_loc) = b_loc(N_loc) + lambda * 0.0_dp  ! u(L,t) = 0
    end if
    
    ! Résolution du système linéaire avec gradient conjugué
    call conjugate_gradient(x_loc, b_loc, a, b, c, N_loc, myrank, nprocs, &
                           left_neighbor, right_neighbor)
    
    ! Mise à jour de la solution
    u_new_loc(1:N_loc) = x_loc(:)
    
    deallocate(b_loc, x_loc)
  end subroutine implicit_step

  ! Algorithme du gradient conjugué pour résoudre le système linéaire
  subroutine conjugate_gradient(x, b, a, b_diag, c, N_loc, myrank, nprocs, &
                               left_neighbor, right_neighbor)
    integer, intent(in) :: N_loc, myrank, nprocs, left_neighbor, right_neighbor
    real(dp), intent(in) :: b(N_loc), a, b_diag, c
    real(dp), intent(out) :: x(N_loc)
    
    real(dp), allocatable :: r(:), p(:), Ap(:)
    real(dp) :: rsold_loc, rsnew_loc, rsold, rsnew
    real(dp) :: pAp_loc, pAp, alpha_cg, beta
    integer :: iter, ierr
    
    ! Allocation des vecteurs temporaires
    allocate(r(N_loc), p(N_loc), Ap(N_loc))
    
    ! Initialisation
    x(:) = b(:)  ! Solution initiale x0 = b
    
    ! Calcul de r0 = b - A*x0
    call compute_Ax(x, Ap, a, b_diag, c, N_loc, myrank, nprocs, left_neighbor, right_neighbor)
    r(:) = b(:) - Ap(:)
    p(:) = r(:)
    
    ! Calculer rsold = <r,r>
    rsold_loc = sum(r(:) * r(:))
    call MPI_ALLREDUCE(rsold_loc, rsold, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! Boucle principale du gradient conjugué
    do iter = 1, max_iter_cg
      ! Calculer A*p
      call compute_Ax(p, Ap, a, b_diag, c, N_loc, myrank, nprocs, left_neighbor, right_neighbor)
      
      ! Calculer pAp = <p, A*p>
      pAp_loc = sum(p(:) * Ap(:))
      call MPI_ALLREDUCE(pAp_loc, pAp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      
      ! Calculer alpha = rsold / pAp
      alpha_cg = rsold / pAp
      
      ! Mettre à jour x et r
      x(:) = x(:) + alpha_cg * p(:)
      r(:) = r(:) - alpha_cg * Ap(:)
      
      ! Calculer rsnew = <r,r>
      rsnew_loc = sum(r(:) * r(:))
      call MPI_ALLREDUCE(rsnew_loc, rsnew, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      
      ! Vérifier la convergence
      if (sqrt(rsnew) < epsilon) exit
      
      ! Calculer beta = rsnew / rsold
      beta = rsnew / rsold
      
      ! Mettre à jour p
      p(:) = r(:) + beta * p(:)
      
      ! Mise à jour pour la prochaine itération
      rsold = rsnew
    end do
    
    deallocate(r, p, Ap)
  end subroutine conjugate_gradient

  ! Calcule le produit matrice-vecteur Ax avec communication aux frontières
  subroutine compute_Ax(x, Ax, a, b, c, N_loc, myrank, nprocs, left_neighbor, right_neighbor)
    integer, intent(in) :: N_loc, myrank, nprocs, left_neighbor, right_neighbor
    real(dp), intent(in) :: x(N_loc), a, b, c
    real(dp), intent(out) :: Ax(N_loc)
    real(dp) :: send_left, send_right, recv_left, recv_right
    integer :: i, ierr
    
    ! Communication pour échanger les valeurs aux frontières
    send_left = x(1)
    send_right = x(N_loc)
    
    call MPI_SENDRECV(send_left, 1, MPI_DOUBLE_PRECISION, left_neighbor, 0, &
                      recv_right, 1, MPI_DOUBLE_PRECISION, right_neighbor, 0, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call MPI_SENDRECV(send_right, 1, MPI_DOUBLE_PRECISION, right_neighbor, 1, &
                      recv_left, 1, MPI_DOUBLE_PRECISION, left_neighbor, 1, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    
    ! Calculer A*x pour les points intérieurs
    do i = 2, N_loc-1
      Ax(i) = b * x(i) + a * x(i-1) + c * x(i+1)
    end do
    
    ! Cas spéciaux pour les points aux frontières du domaine local
    if (myrank == 0) then
      ! Premier processus: premier point avec condition limite à gauche (u(0,t)=0)
      Ax(1) = b * x(1) + c * x(2)
    else
      ! Autres processus: premier point avec valeur reçue de gauche
      Ax(1) = b * x(1) + a * recv_left + c * x(2)
    end if
    
    if (myrank == nprocs-1) then
      ! Dernier processus: dernier point avec condition limite à droite (u(L,t)=0)
      Ax(N_loc) = b * x(N_loc) + a * x(N_loc-1)
    else
      ! Autres processus: dernier point avec valeur reçue de droite
      Ax(N_loc) = b * x(N_loc) + a * x(N_loc-1) + c * recv_right
    end if
  end subroutine compute_Ax

end module solver_mod