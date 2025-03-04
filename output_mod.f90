! filepath: /home/yanel/PA/fortran/first/exercices/heat/output_mod.f90
! -----------------------------------------------------------------------
! Module pour les entrées/sorties et la visualisation
! -----------------------------------------------------------------------

module output_mod
  use mpi
  use parameters_mod, only: dp, N, dx
  implicit none
  private
  
  ! Procédures exportées
  public :: write_solution, print_info
  
contains

  subroutine print_info(step, current_time, myrank)
    integer, intent(in) :: step, myrank
    real(dp), intent(in) :: current_time
    
    ! Afficher des informations de progression (uniquement processus 0)
    if (myrank == 0 .and. mod(step, 10) == 0) then
        print '(A,I5,A,F8.6)', 'Pas de temps ', step, ', t = ', current_time
    end if
  end subroutine print_info
  
  subroutine ensure_directory_exists(myrank)
    integer, intent(in) :: myrank
    logical :: dir_exists
    integer :: ierr
    character(len=100) :: cmd
    
    ! Seulement le processus 0 crée le répertoire
    if (myrank == 0) then
      inquire(file="results/.", exist=dir_exists)
      if (.not. dir_exists) then
        ! Créer le répertoire results s'il n'existe pas
        cmd = "mkdir -p results"
        call execute_command_line(cmd, exitstat=ierr)
        if (ierr /= 0) then
          print *, "Erreur lors de la création du répertoire 'results'"
        end if
      end if
    end if
    
    ! Synchroniser tous les processus
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  end subroutine ensure_directory_exists
  
  subroutine write_solution(u_loc, N_loc, global_start, global_end, step, myrank, nprocs)
    integer, intent(in) :: N_loc, global_start, global_end, step, myrank, nprocs
    real(dp), intent(in) :: u_loc(0:N_loc+1)
    
    integer :: io_unit, i, ierr
    character(len=100) :: filename
    real(dp), allocatable :: u_global(:), x_global(:)
    
    ! S'assurer que le répertoire results existe
    call ensure_directory_exists(myrank)
    
    ! Rassembler la solution complète sur processus 0 pour l'écriture
    if (myrank == 0) then
        allocate(u_global(0:N), x_global(0:N))
        u_global(0) = 0.0_dp  ! Condition limite à gauche
        u_global(N) = 0.0_dp  ! Condition limite à droite
        do i = 0, N
            x_global(i) = i * dx
        end do
    end if
    
    ! Rassembler les données des processus
    call gather_solution(u_loc, N_loc, global_start-1, global_end-1, u_global, N-1, myrank, nprocs)
    
    ! Écrire la solution
    if (myrank == 0) then
        write(filename, '(A,I5.5,A)') 'results/heat_solution_', step, '.dat'
        open(newunit=io_unit, file=trim(filename), status='replace')
        write(io_unit, '(A)') '# x      u(x,t)'
        do i = 0, N
            write(io_unit, '(2F16.8)') x_global(i), u_global(i)
        end do
        close(io_unit)
        deallocate(u_global, x_global)
    end if
  end subroutine write_solution
  
  subroutine gather_solution(u_loc, N_loc, global_start, global_end, u_global, N_interior, myrank, nprocs)
    integer, intent(in) :: N_loc, global_start, global_end, N_interior, myrank, nprocs
    real(dp), intent(in) :: u_loc(0:N_loc+1)
    real(dp), intent(inout) :: u_global(0:N_interior+1)
    integer :: recv_counts(nprocs), displs(nprocs), i, ierr
    real(dp), allocatable :: send_buf(:), recv_buf(:)
    
    ! Chaque processus prépare les données à envoyer (uniquement les points intérieurs)
    allocate(send_buf(N_loc))
    send_buf(:) = u_loc(1:N_loc)
    
    ! Processus 0 collecte les tailles de chaque processus
    call MPI_GATHER(N_loc, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    
    ! Processus 0 calcule les déplacements
    if (myrank == 0) then
        allocate(recv_buf(N_interior))
        displs(1) = 0
        do i = 2, nprocs
            displs(i) = displs(i-1) + recv_counts(i-1)
        end do
    end if
    
    ! Rassemblement des données
    call MPI_GATHERV(send_buf, N_loc, MPI_DOUBLE_PRECISION, &
                    recv_buf, recv_counts, displs, MPI_DOUBLE_PRECISION, &
                    0, MPI_COMM_WORLD, ierr)
    
    ! Processus 0 copie les données dans u_global
    if (myrank == 0) then
        u_global(1:N_interior) = recv_buf(:)
        deallocate(recv_buf)
    end if
    
    deallocate(send_buf)
  end subroutine gather_solution

end module output_mod