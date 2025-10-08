program iar_bench

  use mpi
  use tem_status_module, only: tem_status_type, &
    &                          tem_status_communicate_delayed

  implicit none

  character(len=20) :: nTest_arg
  character(len=40) :: perfFile

  integer, parameter :: rk = selected_real_kind(15)
  integer, parameter :: long_k = selected_int_kind(15)
  integer, parameter :: rk_mpi = MPI_DOUBLE_PRECISION
  integer, parameter :: long_k_mpi = MPI_INTEGER8

  integer, parameter :: totaliter = 10000
  integer, parameter :: samplelen = 1000
  integer, parameter :: test_interval = 1000
  integer :: nTests = 1
  integer :: nRepetitions = totaliter

  integer :: iter
  integer :: intra
  integer :: iError
  integer :: myRank
  integer :: nProcs
  integer :: stat(MPI_STATUS_SIZE)
  integer(kind=long_k) :: res = 0_long_k

  real(kind=rk) :: start_time, stop_time
  real(kind=rk) :: start_comp, stop_comp
  real(kind=rk) :: start_bar, bar_time
  real(kind=rk) :: max_bt
  real(kind=rk) :: compute_time, running_time
  real(kind=rk) :: noncomp_time
  real(kind=rk) :: min_ct, min_rt, min_nc
  real(kind=rk) :: max_ct, max_rt, max_nc
  real(kind=rk) :: avg_ct, avg_rt, avg_nc

  logical :: iar_completed = .false.

  type(tem_status_type) :: health

  call MPI_Init(iError)

  call MPI_Comm_size(MPI_COMM_WORLD, nProcs, iError)
  call MPI_Comm_rank(MPI_COMM_WORLD, myRank, iError)

  if (command_argument_count() > 0) then
    call get_command_argument(1, nTest_arg)
    read(nTest_arg, *) nTests
    if (mod(totaliter, nTests) /= 0) then
       nTests = totaliter / (totaliter/nTests)
    end if
    nRepetitions = totaliter/nTests
  end if
  compute_time = 0.0_rk

  start_bar = MPI_Wtime()
  call MPI_Barrier(MPI_COMM_WORLD, iError)
  start_time = MPI_Wtime()
  do iter = 1,nRepetitions
    iar_completed = .false.

    start_comp = MPI_Wtime()
    do intra=1,nTests
      call compute(test_interval, res)
      if (.not. iar_completed) then
        call MPI_Test(health%check_request, iar_completed, stat, iError)
      end if
    end do
    health%bits(1) = res > 1000
    stop_comp = MPI_Wtime()
    compute_time = compute_time + stop_comp - start_comp
    call tem_status_communicate_delayed( me   = health,        &
      &                                  comm = MPI_COMM_WORLD )
  end do
  stop_time = MPI_Wtime()

  running_time = stop_time - start_time
  noncomp_time = running_time - compute_time
  bar_time = start_time - start_bar

  call MPI_Reduce(compute_time, min_ct, 1, rk_mpi, MPI_MIN, 0, MPI_COMM_WORLD, iError)
  call MPI_Reduce(compute_time, max_ct, 1, rk_mpi, MPI_MAX, 0, MPI_COMM_WORLD, iError)
  call MPI_Reduce(compute_time, avg_ct, 1, rk_mpi, MPI_SUM, 0, MPI_COMM_WORLD, iError)

  call MPI_Reduce(running_time, min_rt, 1, rk_mpi, MPI_MIN, 0, MPI_COMM_WORLD, iError)
  call MPI_Reduce(running_time, max_rt, 1, rk_mpi, MPI_MAX, 0, MPI_COMM_WORLD, iError)
  call MPI_Reduce(running_time, avg_rt, 1, rk_mpi, MPI_SUM, 0, MPI_COMM_WORLD, iError)

  call MPI_Reduce(noncomp_time, min_nc, 1, rk_mpi, MPI_MIN, 0, MPI_COMM_WORLD, iError)
  call MPI_Reduce(noncomp_time, max_nc, 1, rk_mpi, MPI_MAX, 0, MPI_COMM_WORLD, iError)
  call MPI_Reduce(noncomp_time, avg_nc, 1, rk_mpi, MPI_SUM, 0, MPI_COMM_WORLD, iError)

  call MPI_Reduce(bar_time, max_bt, 1, rk_mpi, MPI_MAX, 0, MPI_COMM_WORLD, iError)

  if (myRank == 0) then
    avg_ct = avg_ct / nProcs
    avg_rt = avg_rt / nProcs
    avg_nc = avg_nc / nProcs
    write(perfFile, '(A7,I7.7)') 'ncavg.p', nProcs
    write(*,*) health%bits
    write(*,*) ''
    write(*,*) 'Max. Barrier time:', max_bt
    write(*,*) 'Running time: min=', min_rt, '; max=', max_rt, '; avg=', avg_rt
    write(*,*) 'Compute time: min=', min_ct, '; max=', max_ct, '; avg=', avg_ct
    write(*,*) 'NonComp time: min=', min_nc/nRepetitions, &
      &                    '; max=', max_nc/nRepetitions, &
      &                    '; avg=', avg_nc/nRepetitions
    open(22, FILE=perfFile, POSITION='APPEND', FORM='FORMATTED', ACTION='WRITE')
    write(22,*) nTests, avg_nc/nRepetitions
    close(22)
  end if
  call MPI_Finalize(iError)

contains

  subroutine compute(nComputations, res)
    integer, intent(in) :: nComputations
    integer(kind=long_k), intent(out) :: res


    integer :: iComp
    real(kind=rk) :: x(samplelen), y(samplelen)
    real(kind=rk) :: r(samplelen)

    res = 0_long_k

    do iComp=1,nComputations
      call random_number(x)
      call random_number(y)
      r = sqrt(x**2 + y**2)
      res = res + count(r <= 1.0_rk)
    end do
  end subroutine compute

end program iar_bench
