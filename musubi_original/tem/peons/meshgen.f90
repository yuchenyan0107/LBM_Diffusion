PROGRAM main 
  use mpi
  integer, parameter                   :: long_k = selected_int_kind(15)
  integer(kind=long_k)                 :: Elem, i, elemOffset
  integer(kind=long_k), allocatable    :: buf(:) 
  integer(kind=MPI_OFFSET_KIND)        :: filesize, disp
  integer                              :: ierr, myrank, mysize, fh, etype, ftype, remainder, pos_Elem
  integer                              :: maxElem_write, write_Elem, long_k_mpi, loc_Elem, nChunk, comm
  CHARACTER(len=32)                    :: arg
  long_k_mpi = MPI_INTEGER8

  ! Generates a mesh for MPIIO testcase
  ! number of Elements is 700M on default, if not given as Argument
   

  ! MPI init
  call MPI_INIT(ierr) 
  comm=MPI_COMM_WORLD
  call MPI_COMM_RANK(comm, myrank, ierr) 
  call MPI_COMM_SIZE(comm, mysize, ierr) 

  ! Check Arguments for number of elements to create mesh with
  if (myrank==0)then
    write(*,*) "Generating Testmesh for MPIIO test program."
    write(*,*) "------------------------------------------------------------------------------------"
  endif
  ! Check if Number of Elements are given as argument, default = 720M
  if ( COMMAND_ARGUMENT_COUNT() == 1 ) then
    CALL get_command_argument(1, arg)
    read(arg(:),' (I16)') Elem 
    if (myrank==0)then
      write(*,*) "Generate mesh with", Elem ,' Elements.'
    endif
  else
    ! if no argument given, use 720M Elements as default
    Elem = 720000000 
    if (myrank==0)then
      write(*,*) "No Argument given use default of " ,Elem,' Elements!'
      write(*,*) "To generate a mesh with a different amount of Elements pass the amount as Argument."
    endif
  end if

  ! calculate local Element count
  loc_Elem= Elem/mysize

  ! Check if local Elements are below max value 
  if(loc_Elem*16_long_k > huge(loc_Elem))then
    if (myrank==0)then
      write(*,*) "Exceed maximum of", huge(loc_Elem)/16, "Elements per process, program terminatated!"
      call MPI_abort(comm,1,ierr)
    endif
  endif

  ! Get amount of neccessery chunks cause of MPI 2GB write limit per write
  maxElem_write=huge((loc_Elem)/16)-1
  nChunk = int(loc_Elem/maxElem_write)
  if (mod(loc_Elem, maxElem_write) == 0 ) then
    nChunk=nChunk-1
  endif

  ! Distribute remaining Elements
  remainder = int(mod(Elem, mySize))
  elemOffset = loc_Elem * myRank + min(myRank, remainder)
  if (myrank < mod( Elem, mysize) ) then      
    loc_Elem = loc_Elem + 1
  endif 

  !allocate buffer
  allocate(buf(loc_Elem*2))

  ! fill in some values
  do i = 0, loc_Elem-1
    if (myrank==0) then
      buf(i*2+1) = i 
      buf(i*2+2) = i 
    else
      buf(i*2+1) = elemOffset + i 
      buf(i*2+2) = elemOffset + i 
    endif
  enddo

  ! Create a contiguous type to describe the vector per element
  call MPI_TYPE_CONTIGUOUS( 2, long_k_mpi, etype, iError )
  call MPI_TYPE_COMMIT(etype, iError )      
  call check_mpi_error( ierr, 'etype', comm)
  call MPI_TYPE_CONTIGUOUS( loc_Elem, etype, ftype, iError )
  call MPI_TYPE_COMMIT(ftype, iError )      
  call check_mpi_error( ierr, 'ftype', comm)

  ! Calculate Displacement                   
  disp = elemOffset * 16_MPI_OFFSET_KIND 

  ! MPI IO 
  call MPI_FILE_OPEN(comm, 'mesh.out', MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, ierr) 
  call check_mpi_error( ierr, 'Open', comm)
  call MPI_FILE_SET_VIEW(fh, disp , etype, ftype, 'native', MPI_INFO_NULL, ierr)
  call check_mpi_error( ierr, 'Setview', comm)
  ! Dump several time cause of 2GB MPI limit per write call
  do i = 0, nChunk
    write_Elem = min((loc_Elem-i*maxElem_write), maxElem_write)
    pos_Elem = 2*i*maxElem_write+1 
    call MPI_FILE_WRITE_ALL(fh, buf(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, ierr) 
  enddo
  call check_mpi_error( ierr , 'Write', comm)
  call MPI_BARRIER(comm,IERR)
  call MPI_FILE_CLOSE(fh, ierr) 
  call check_mpi_error( ierr , 'Close', comm)
  if (myrank==0)then
    write(*,*) "Successfully generated mesh.out"
  endif
  call MPI_FINALIZE(ierr) 
   
contains
 ! MPI ERROR CHECK
 subroutine check_mpi_error( iError, event_string, comm )
   integer, intent(in) :: iError, comm
   character(len=*), intent(in) :: event_string
   character(len=100) :: IOError
   integer :: resultlen = 100
   integer :: ErrErr
 
   if (iError /= MPI_SUCCESS) then
     call MPI_ERROR_STRING(iError, IOError, resultlen, ErrErr)
     write(*,*) 'MPI Error when '//trim(event_string),': ' &
       &                 //trim(IOError)
     call mpi_abort(comm,1,ErrErr)
   end if
 end subroutine check_mpi_error
END PROGRAM main 
