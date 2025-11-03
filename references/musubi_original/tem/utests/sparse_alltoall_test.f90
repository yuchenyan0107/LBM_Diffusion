! Copyright (c) 2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
!mpi!nprocs = 4
program sparse_alltoall_test
  use mpi
  use tem_sparse_comm_module, only: tem_sparse_alltoall_int

  implicit none

  ! Assume we have 4 procs: 0 - 1 - 2 - 3
  ! Each proc send data to their direct neighbors
  ! the data is a single integer whose value is the targeting proc, i.e.
  ! proc 0 sends value 1 to proc 1
  ! proc 1 sends value 0 to proc 0, value 2 to proc 2
  ! proc 2 sends value 1 to proc 1, value 3 to proc 3
  ! proc 3 sends value 2 to proc 2
  ! After communication, each proc should have the following data:
  ! proc 0: 0
  ! proc 1: 1 1
  ! proc 2: 2 2
  ! proc 3: 3

  integer :: iError
  integer :: myrank, nProcs
  logical :: local_success, success
  ! integer :: rstat(MPI_STATUS_SIZE)

  integer, allocatable :: targets(:)
  integer, allocatable :: send_buffer(:)
  integer, allocatable :: sources(:)
  integer, allocatable :: recv_buffer(:)

  call MPI_Init(iError)
  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, iError)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, IERROR)

  if ( myrank == 0 ) then
    ! has one neighbor
    allocate( targets(1) )
    allocate( send_buffer(1) )
    targets(1)     = myrank + 1
    send_buffer(1) = myrank + 1
  else if ( myrank == (nProcs-1) ) then
    ! has one neighbor
    allocate( targets(1) )
    allocate( send_buffer(1) )
    targets(1)     = myrank - 1
    send_buffer(1) = myrank - 1
  else
    ! has two neighbor
    allocate( targets(2) )
    allocate( send_buffer(2) )
    targets(1)     = myrank - 1
    targets(2)     = myrank + 1
    send_buffer(1) = myrank - 1
    send_buffer(2) = myrank + 1
  end if

  call tem_sparse_alltoall_int( targets, send_buffer, &
    &                           sources, recv_buffer, MPI_COMM_WORLD)

  if ( myrank == 0 .or. myrank == (nProcs-1) ) then
    local_success = (recv_buffer(1) == myrank)
  else
    local_success = ((recv_buffer(1) == myrank) .and. (recv_buffer(2) == myrank))
  end if

  if ( .not. local_success ) then
    write(*,*) "Rank: ", myrank, " sources: ", sources
    write(*,*) "Rank: ", myrank, " recv_bu: ", recv_buffer
  end if

  call MPI_REDUCE(local_success, success, 1, MPI_LOGICAL, MPI_LAND, 0, MPI_COMM_world, IERROR)

  if ( myrank == 0 ) then
    if ( success ) then
      write(*,*) "PASSED"
    else
      write(*,*) "FAILED"
    end if
  end if

  call MPI_Finalize(iError)

end program sparse_alltoall_test
