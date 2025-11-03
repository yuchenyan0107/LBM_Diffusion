program checkTest
  use mpi
  use distCheck_module, only: checksum

  implicit none

  real :: somedat(100)
  real, allocatable :: alldat(:,:)
  integer :: iError, myRank, rl
  integer :: nProcs

  call MPI_Init(iError)

  call MPI_Comm_Rank(MPI_COMM_WORLD, myRank, iError)
  call MPI_Comm_Size(MPI_COMM_WORLD, nProcs, iError)
  allocate(alldat(100, nProcs))

  somedat = myRank

  write(*,'(a)') checksum(somedat, shape(somedat), MPI_COMM_WORLD)

  call MPI_Gather(somedat, 100, MPI_Real, alldat, 100, MPI_Real, 0, &
    &             MPI_COMM_WORLD, iError)
  if (myRank == 0) then
    write(*,*) 'Gathered on Proc 0:', checksum(alldat, shape(alldat), &
      &                                        MPI_COMM_SELF)
  end if

  inquire(iolength=rl) somedat
  if (myRank == 0) then
    open(file = "ref.dat", action = "write", status = "replace", &
         form = "unformatted", access = "direct", unit = 23, recl=rl)
    write(23,rec=myRank+1) somedat
    call MPI_Barrier(MPI_COMM_WORLD, iError)
  else
    call MPI_Barrier(MPI_COMM_WORLD, iError)
    open(file = "ref.dat", action = "write", status = "old", &
         form = "unformatted", access = "direct", unit = 23, recl=rl)
    write(23,rec=myRank+1) somedat
  end if
  close(23)

  call MPI_Finalize(iError)

end program checkTest
