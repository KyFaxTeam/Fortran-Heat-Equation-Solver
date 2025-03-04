! -----------------------------------------------------------------------
! Programme principal: Résolution de l'équation de la chaleur en 1D avec MPI
! -----------------------------------------------------------------------
! Ce programme utilise plusieurs modules pour résoudre l'équation de la 
! chaleur en une dimension spatiale avec une méthode implicite et MPI.
! -----------------------------------------------------------------------

program heat_solver
  use mpi
  use parameters_mod, only: dp, N, dt, t_final, alpha, dx, lambda
  use domain_mod, only: setup_domain, compute_local_indices, exchange_boundaries
  use solver_mod, only: implicit_step, conjugate_gradient
  use output_mod, only: write_solution, print_info
  implicit none

  ! Variables MPI
  integer :: ierr, myrank, nprocs
  
  ! Variables pour le domaine local
  integer :: N_loc, global_start, global_end
  integer :: left_neighbor, right_neighbor
  
  ! Tableaux locaux pour la solution
  real(dp), allocatable :: u_loc(:), u_new_loc(:)
  
  ! Variables pour la boucle temporelle
  integer :: step, max_steps
  real(dp) :: current_time
  
  ! Initialisation MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  
  ! Configuration du domaine et calcul des indices locaux
  call compute_local_indices(N, nprocs, myrank, N_loc, global_start, global_end, &
                            left_neighbor, right_neighbor)
  
  ! Allouer les tableaux locaux
  allocate(u_loc(0:N_loc+1), u_new_loc(0:N_loc+1))
  
  ! Initialiser le domaine avec la condition initiale
  call setup_domain(u_loc, N_loc, global_start, global_end)
  
  ! Paramètres temporels
  current_time = 0.0_dp
  max_steps = nint(t_final / dt)
  
  ! Boucle temporelle principale
  do step = 1, max_steps
      ! Mettre à jour le temps courant
      current_time = current_time + dt
      
      ! Effectuer un pas de temps implicite avec gradient conjugué
      call implicit_step(u_loc, u_new_loc, N_loc, alpha, dt, myrank, nprocs, &
                         left_neighbor, right_neighbor)
      
      ! Échange des valeurs aux frontières
      call exchange_boundaries(u_new_loc, N_loc, myrank, nprocs, &
                              left_neighbor, right_neighbor)
      
      ! Mettre à jour u_loc pour le pas de temps suivant
      u_loc(:) = u_new_loc(:)
      
      ! Afficher des informations de progression
      call print_info(step, current_time, myrank)
      
      ! Écrire la solution à certains pas de temps
      if (mod(step, 10) == 0) then
          call write_solution(u_loc, N_loc, global_start, global_end, step, myrank, nprocs)
      end if
  end do
  
  ! Libération de la mémoire
  deallocate(u_loc, u_new_loc)
  
  ! Finalisation MPI
  call MPI_FINALIZE(ierr)

end program heat_solver