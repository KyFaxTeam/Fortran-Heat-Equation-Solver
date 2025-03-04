! -----------------------------------------------------------------------
! Module pour la gestion du domaine et la décomposition parallèle
! -----------------------------------------------------------------------

module domain_mod
  use mpi
  use parameters_mod, only: dp, N, L, dx, pi
  implicit none
  private
  
  ! Procédures exportées
  public :: compute_local_indices, setup_domain, exchange_boundaries
  
contains

  ! Calcule les indices locaux pour la décomposition du domaine
  subroutine compute_local_indices(N_global, nprocs, myrank, N_loc, global_start, global_end, &
                                  left_neighbor, right_neighbor)
    integer, intent(in) :: N_global, nprocs, myrank
    integer, intent(out) :: N_loc, global_start, global_end
    integer, intent(out) :: left_neighbor, right_neighbor
    integer :: i, base_size, remainder
    
    ! Calcul de la taille locale (distribution équilibrée)
    base_size = (N_global-1) / nprocs
    remainder = mod(N_global-1, nprocs)
    
    ! Si le rang est inférieur au reste, ce processus reçoit un point de plus
    if (myrank < remainder) then
      N_loc = base_size + 1
    else
      N_loc = base_size
    end if
    
    ! Calcul des indices globaux pour ce processus
    global_start = 1  ! Point de départ (indice global)
    do i = 0, myrank-1
      global_start = global_start + base_size
      if (i < remainder) global_start = global_start + 1
    end do
    global_end = global_start + N_loc - 1
    
    ! Identification des voisins
    left_neighbor = myrank - 1
    right_neighbor = myrank + 1
    if (myrank == 0) left_neighbor = MPI_PROC_NULL
    if (myrank == nprocs-1) right_neighbor = MPI_PROC_NULL
  end subroutine compute_local_indices

  ! Initialise le domaine local avec la condition initiale
  subroutine setup_domain(u_loc, N_loc, global_start, global_end)
    integer, intent(in) :: N_loc, global_start, global_end
    real(dp), intent(out) :: u_loc(0:N_loc+1)
    integer :: i
    real(dp) :: x
    
    ! Initialisation avec la condition initiale u(x,0) = sin(πx)
    do i = 1, N_loc
      x = (global_start + i - 1) * dx
      u_loc(i) = sin(pi * x)
    end do
    
    ! Initialiser les points fantômes (seront mis à jour par communication)
    u_loc(0) = 0.0_dp
    u_loc(N_loc+1) = 0.0_dp
  end subroutine setup_domain

  ! Échange les valeurs aux frontières entre les processus
  subroutine exchange_boundaries(u_loc, N_loc, myrank, nprocs, left_neighbor, right_neighbor)
    integer, intent(in) :: N_loc, myrank, nprocs, left_neighbor, right_neighbor
    real(dp), intent(inout) :: u_loc(0:N_loc+1)
    real(dp) :: send_left, send_right, recv_left, recv_right
    integer :: ierr
    
    ! Valeurs à envoyer aux voisins
    send_left = u_loc(1)
    send_right = u_loc(N_loc)
    
    ! Échange avec les voisins (non-bloquant pour plus d'efficacité)
    call MPI_SENDRECV(send_left, 1, MPI_DOUBLE_PRECISION, left_neighbor, 0, &
                      recv_right, 1, MPI_DOUBLE_PRECISION, right_neighbor, 0, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call MPI_SENDRECV(send_right, 1, MPI_DOUBLE_PRECISION, right_neighbor, 1, &
                      recv_left, 1, MPI_DOUBLE_PRECISION, left_neighbor, 1, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    
    ! Mise à jour des points fantômes
    if (right_neighbor /= MPI_PROC_NULL) u_loc(N_loc+1) = recv_right
    if (left_neighbor /= MPI_PROC_NULL) u_loc(0) = recv_left
    
    ! Appliquer les conditions aux limites de Dirichlet
    if (myrank == 0) u_loc(0) = 0.0_dp  ! u(0,t) = 0
    if (myrank == nprocs-1) u_loc(N_loc+1) = 0.0_dp  ! u(L,t) = 0
  end subroutine exchange_boundaries

end module domain_mod