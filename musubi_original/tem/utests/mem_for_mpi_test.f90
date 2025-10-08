! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this
! list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
! this list of conditions and the following disclaimer in the documentation
! and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
! DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!mpi!nprocs = 2 
program mem_for_mpi_test
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  use hvs_sizeof_module, only: c_sizeof
  use mpi
  use mem_for_mpi_module

  implicit none

  integer          :: iError
  type(c_ptr)      :: buffer
  integer, pointer :: fordat(:)
  integer          :: intsize
  integer          :: myrank
  integer          :: rstat(MPI_STATUS_SIZE)
  logical          :: success

  call MPI_Init(iError)

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, iError)

  if (myrank == 0) then
    write(*,*) 'Starting test for alloc_mpif_mem'
  end if
  intsize = c_sizeof(iError)

  ! Allocate 10 integers:
  call alloc_mpif_mem( asize   = intsize*10_MPI_ADDRESS_KIND, &
    &                  baseptr = buffer,                      &
    &                  ierr    = iError                       )
  if (myrank == 0) then
    write(*,*) 'Allocated memory with iError=', iError

    write(*,*) 'Converting C pointer to Fortran array'
  end if
  call c_f_pointer(buffer, fordat, [10])

  if (myrank == 0) then
    write(*,*) 'Assigning a value to the buffer with size:', size(fordat)
    fordat = 42
    write(*,*) 'Sending data from rank 0 to rank 1'
    call MPI_Send(fordat, 10, MPI_Integer, 1, 23, MPI_COMM_WORLD, iError)
    call MPI_Recv(success, 1, MPI_Logical, 1, 32, MPI_COMM_WORLD, rstat, iError)
    if (success) then
      write(*,*) 'PASSED'
    else
      write(*,*) 'FAILED'
    end if
  else
    call MPI_Recv(fordat, 10, MPI_Integer, 0, 23, MPI_COMM_WORLD, rstat, iError)
    if (all(fordat == 42)) then
      nullify(fordat)
      call free_mpif_mem(buffer)
      success = .true.
    else
      success = .false.
    end if
    call MPI_Send(success, 1, MPI_Logical, 0, 32, MPI_COMM_WORLD, iError)
  end if

  call MPI_Finalize(iError)

end program mem_for_mpi_test
