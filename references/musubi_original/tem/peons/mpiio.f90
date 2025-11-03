program MPIIO
      use MPI
      implicit none
      integer, parameter                :: long_k = selected_int_kind(15)
      integer, parameter                :: rk = selected_real_kind(15)
      CHARACTER(len=30)                 :: arg, filename, outfile
      character(len=30)                 :: mpicase, timing_filename
      real(kind=rk)                     :: start
      real(kind=rk)                     :: glob_iotime
      integer(kind=rk)                  :: nElems, elemOffset
      integer(kind=long_k), allocatable :: buffer(:)
      integer(kind=MPI_OFFSET_KIND)     :: displacement, filesize
      integer                           :: nlocElems, nParts, comm, rank
      integer                           :: fh, fho, etype, ftype, method
      integer                           :: iError, ierr, nChunk, i
      integer                           :: maxElem_write, write_Elem, pos_Elem 

      ! write information about this program
      if (rank==0)then
        write(*,*) "This test program tries read and write treelm mesh file"
      endif

      ! Name of files
      filename         = 'mesh.out'
      outfile          = 'outmesh.out'
      timing_filename  = 'timing.res'

      ! Init MPI
      call MPI_Init(iError)
      comm             = MPI_COMM_WORLD
      Call MPI_COMM_SIZE(comm,nParts,iError)
      Call MPI_COMM_rank(comm,rank,iError)

      ! Check arguments for case
      if ( COMMAND_ARGUMENT_COUNT() == 1 ) then
        CALL get_command_argument(1, arg)
        IF ((arg == '1').or. (arg == '2').or.(arg == '3').or.(arg == '4')) then 
          read(arg(1:1),' (I2)') method
        else
          if (rank==0)then
            write(*,*) "Please select a method 1 - 4"
            write(*,*) "Use 1 to test MPI_FILE_READ        / MPI_FILE_WRITE"
            write(*,*) "Use 2 to test MPI_FILE_READ_AT     / MPI_FILE_WRITE_AT "
            write(*,*) "Use 3 to test MPI_FILE_READ_ALL    / MPI_FILE_WRITE_ALL"
            write(*,*) "Use 4 to test MPI_FILE_READ_AT_ALL / MPI_FILE_WRITE_AT_ALL"
            call mpi_abort(comm,1,ierr)
          endif
        endif
      else
        if (rank==0)then
          write(*,*) "Please pass the argument for choosen method:"
          write(*,*) "Available methods 1 - 4:"
          write(*,*) "Use 1 to test MPI_FILE_READ        / MPI_FILE_WRITE"
          write(*,*) "Use 2 to test MPI_FILE_READ_AT     / MPI_FILE_WRITE_AT "
          write(*,*) "Use 3 to test MPI_FILE_READ_ALL    / MPI_FILE_WRITE_ALL"
          write(*,*) "Use 4 to test MPI_FILE_READ_AT_ALL / MPI_FILE_WRITE_AT_ALL"
          call mpi_abort(comm,1,ierr)
        endif
      end if
      
      ! Get neccessery values
      call init_values(nParts,rank,comm,nElems,nlocElems,elemOffset,filename,maxElem_write,nChunk)

      select case ( method )
!---------------------------------------------------------------------------------------------------------------
         case( 1 ) ! Use read/write
           mpicase='read'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_READ"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fh, filename, outfile, comm, ElemOffset, 0)
           call MPI_FILE_SET_VIEW( fh, displacement , etype, ftype , "native", MPI_INFO_NULL, iError )
           call check_mpi_error( iError,"SetView" ) 
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             call MPI_FILE_READ( fh, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, iError )
           enddo
           call check_mpi_error( iError,"Read" ) 
           call MPI_FILE_CLOSE(fh, iError)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)
           
           mpicase='write'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_WRITE"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fho, filename, outfile, comm, ElemOffset, 1)
           call MPI_FILE_SET_VIEW(fho, displacement , etype, ftype, 'native', MPI_INFO_NULL, ierr)
           call check_mpi_error( iError,"SetView" ) 
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             call MPI_FILE_WRITE(fho, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, ierr)
           enddo
           call check_mpi_error( iError,"Write" ) 
           call MPI_FILE_CLOSE(fho, ierr)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)
!---------------------------------------------------------------------------------------------------------------
         case( 2 ) ! Use read_at/write_at
           mpicase='read_at'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_READ_AT"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fh, filename, outfile, comm, ElemOffset, 0)
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             displacement = (ElemOffset+i*maxElem_write)*16_MPI_OFFSET_KIND
             call MPI_FILE_READ_AT(fh, displacement, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, iError )
           enddo
           call MPI_FILE_CLOSE(fh, iError)
           call Finish_timer(start, nParts ,nElems, mpicase, timing_filename, comm, glob_iotime)

           mpicase='write_at'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_WRITE_AT"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fho, filename, outfile, comm, ElemOffset, 1)
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             displacement = (ElemOffset+i*maxElem_write)*16_MPI_OFFSET_KIND
             call MPI_FILE_WRITE_AT(fho, displacement, buffer(pos_Elem:),write_Elem , etype, MPI_STATUS_IGNORE, ierr)
           enddo
           call check_mpi_error( ierr , 'write')
           call MPI_FILE_CLOSE(fho, ierr)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)
!---------------------------------------------------------------------------------------------------------------
        case( 3 ) ! Use read_all/write_all
           mpicase='read_all'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_READ_ALL"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fh, filename, outfile, comm, ElemOffset, 0)
           call MPI_FILE_SET_VIEW( fh, displacement , etype, ftype , "native", MPI_INFO_NULL, iError )
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             call MPI_FILE_READ_ALL( fh, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, iError )
           enddo
           call MPI_FILE_CLOSE(fh, iError)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)

           mpicase='write_all'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_WRITE_ALL"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fho, filename, outfile, comm, ElemOffset, 1)
           call MPI_FILE_SET_VIEW(fho, displacement , etype, ftype, 'native', MPI_INFO_NULL, ierr)
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             call MPI_FILE_WRITE_ALL(fho, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, ierr)
           enddo
           call MPI_FILE_CLOSE(fho, ierr)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)
!---------------------------------------------------------------------------------------------------------------
         case( 4 ) ! Use read_at_all/write_at_all
           mpicase='read_at_all'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_READ_AT_ALL"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fh, filename, outfile, comm, ElemOffset, 0)
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             displacement = (ElemOffset+i*maxElem_write)*16_MPI_OFFSET_KIND
             call MPI_FILE_READ_AT_ALL( fh, displacement, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, iError )
           enddo
           call MPI_FILE_CLOSE(fh, iError)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)
  
           mpicase='write_at_all'
           if ( rank==0 ) write(*,*) "Using MPI_FILE_WRITE_AT_ALL"
           call Start_timer(comm,start)
           call init_MPI(etype, ftype, nlocElems, displacement, fho, filename, outfile, comm, ElemOffset, 1)
           do i = 0, nChunk
             write_Elem = min((nlocElems-i*maxElem_write), maxElem_write)
             pos_Elem = 2*i*maxElem_write+1 
             displacement = (ElemOffset+i*maxElem_write)*16_MPI_OFFSET_KIND
             call MPI_FILE_WRITE_AT_ALL(fho, displacement, buffer(pos_Elem:), write_Elem, etype, MPI_STATUS_IGNORE, ierr)
           enddo
           call MPI_FILE_CLOSE(fho, ierr)
           call Finish_timer(start, nParts ,nElemS, mpicase, timing_filename, comm, glob_iotime)
!---------------------------------------------------------------------------------------------------------------
      end select

  call MPI_Finalize(iError)


contains
  !initiation routines
  subroutine init_values(nParts,myRank,comm,nElems,nlocElems,elemOffset,filename,maxElem_write,nChunk)
    integer, parameter                :: rk = selected_real_kind(15)
    integer(kind=rk),intent(out)      :: nElems, elemOffset 
    integer,intent(in)                :: nParts, comm, myRank
    integer,intent(out)               :: nlocElems, nChunk
    integer                           :: iError, fh, remainder, maxElem_write
    integer(kind=MPI_OFFSET_KIND)     :: filesize
    Character(len=30),intent(in)      :: filename

    ! Open the binary file for MPI I/O (Write)
    call MPI_FILE_OPEN( comm, filename, MPI_MODE_RDONLY,   &
      &                 MPI_INFO_NULL, fh, iError )
    !Get number of Elements
    call MPI_File_get_size(fh, filesize, iError)
    call MPI_FILE_CLOSE(fh, ierr)
    nElems = filesize/16_MPI_OFFSET_KIND
    !! scatter elements equal to all process 
    nlocElems  = nElems/nParts
    remainder = int(mod(nElems, nParts))
    elemOffset = nlocElems * myRank + min(myRank, remainder)
    !! scatter restelements 
    if (rank < mod( nElems, nParts) ) then      
             nlocElems = nlocElems + 1
    endif 
    ! allocationts
    allocate(buffer(nElems*2))

    ! Get amount of neccessery chunks cause of MPI 2GB write limit per write
    maxElem_write=(1.9*1024*1024*1024)/16
    if(nlocElems*16_long_k > huge(nlocElems))then
      write(*,*) "Exceed maximum of 134217727 Elements per process, program terminatated!"
      write(*,*) "Please restart with more processes!"
      call MPI_abort(comm,1,ierr)
    endif
    nChunk = int(nlocElems/maxElem_write)
    if (mod(nlocElems, maxElem_write) == 0 ) then
      nChunk=nChunk-1
    endif
  endsubroutine init_values


  ! Initiate MPIIO, open files and initiate types
  subroutine init_MPI(etype, ftype, nlocElems, displacement, fh, filename, outfile, comm, ElemOffset, icase)
    integer, parameter                :: rk = selected_real_kind(15)
    integer                           :: iError, icase
    integer,intent(in)                :: nlocElems, comm
    integer(kind=rk),intent(in)       :: elemOffset
    integer,intent(out)               :: fh, etype, ftype
    integer(kind=MPI_OFFSET_KIND)     :: displacement
    Character(len=30),intent(in)      :: filename, outfile
    ! Open File for read/write
    if ( icase == 1 ) then
      call MPI_FILE_OPEN( comm, outfile, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, fh, iError )
      call check_mpi_error( ierr , 'open')
    else
      call MPI_FILE_OPEN( comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, fh, iError )
      call check_mpi_error( ierr , 'open')
    endif
    ! Create a contiguous type to describe the vector per element
    call MPI_TYPE_CONTIGUOUS( 2, MPI_INTEGER8, etype, iError )
    call MPI_TYPE_COMMIT(etype, iError )
    ! Create a contiguous type to describe the vector per rank
    call MPI_TYPE_CONTIGUOUS( nlocElems, etype, ftype, iError )
    call MPI_TYPE_COMMIT(ftype, iError )
    ! Calculate Displacement
    displacement = ElemOffset*16_MPI_OFFSET_KIND
  endsubroutine init_MPI


  !Timer controle
  subroutine start_timer(comm, start)
    integer, parameter                :: rk = selected_real_kind(15)
    integer,intent(in)                :: comm
    real(kind=rk),intent(out)         :: start
    call MPI_Barrier(comm, iError)
    start = MPI_Wtime()         
  endsubroutine start_timer


  ! Time measurement
  subroutine Finish_timer(start,nParts, nElems, mpicase, timing_filename, comm, glob_iotime)
    integer, parameter                :: rk = selected_real_kind(15)
    integer,intent(in)                :: comm, nParts 
    integer(kind=rk),intent(in)       :: nElems
    real(kind=rk),intent(in)          :: start
    real(kind=rk)                     :: finish, t_difference  
    real(kind=rk),intent(inout)       :: glob_iotime
    character(30),intent(in)          :: mpicase, timing_filename
    finish       = MPI_Wtime()
    t_difference = finish - start
    call MPI_Reduce(t_difference, glob_iotime, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, iError)
    ! write timings to File
    if (rank == 0) then
      call write_timing( nParts, nElems,  mpicase, glob_iotime , timing_filename)
    end if
  endsubroutine Finish_timer


  !Timing dump
  subroutine write_timing( nParts, nElems, mpicase, glob_iotime , timing_filename)
    ! Dumps timing to timing file 
    integer, parameter                :: rk = selected_real_kind(15)
    integer(kind=rk),intent(in)       :: nElems
    integer                           :: fileunit, stat=0, nParts
    logical                           :: file_exists
    character(30),intent(in)          :: mpicase
    character(30)                     :: timing_filename
    character(1024)                   :: header, output
    real(kind=rk),intent(in)          :: glob_iotime
    header= ''
    output= ''
    write(header,'(a)') '# nProcs'
    write(output,'(I8.0)') nParts
    write(header,'(a,a12)') trim(header),'MPI_case'
    write(output,'(a,1x,a11)') trim(output), trim(mpicase)
    write(header,'(a,1x,a10)') trim(header),'nElems'
    write(output,'(a,1x,I10)') trim(output), nElems
    write(header,'(a,1x,a15)') trim(header),'Time_needed'
    write(output,'(a,1x,E15.7)') trim(output), glob_iotime
    inquire(file=timing_filename,exist=file_exists)
    fileunit=newunit()
    open(unit=fileunit,file=trim(timing_filename),position='append',iostat=stat)
    if( .not. file_exists ) then
      write(fileunit,*) trim(header)
    end if
    write(fileunit,*) trim(output)
    close(fileunit)
  end subroutine write_timing


  function newunit() result(nu)
    integer :: nu
    integer, save :: nu_start = 22
    logical :: connected
  
    nu = nu_start
    inquire(unit=nu, opened=connected)
    do while(connected)
      nu = nu + 1
      inquire(unit=nu, opened=connected)
    end do
    nu_start = nu+1
  end function newunit

 ! MPI ERROR CHECK
 subroutine check_mpi_error( iError, event_string )
   integer, intent(in) :: iError
   character(len=*), intent(in) :: event_string
   character(len=100) :: IOError
   integer :: resultlen = 100
   integer :: ErrErr
 
   if (iError /= MPI_SUCCESS) then
     call MPI_ERROR_STRING(iError, IOError, resultlen, ErrErr)
     write(*,*) 'MPI Error when '//trim(event_string),': ' &
       &                 //trim(IOError)
   end if
 end subroutine check_mpi_error

end program MPIIO 
