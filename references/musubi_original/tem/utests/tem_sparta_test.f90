! Copyright (c) 2016 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016,2019,2022 Harald Klimach <harald.klimach@dlr.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
!mpi!nprocs = 5
program tem_sparta_test
  use mpi
  use env_module,               only: rk, init_env, fin_env, long_k
  use tem_utestEnv_module,  only: cubeconf
  use tem_sparta_module, only: tem_sparta_type, tem_balance_sparta, &
    &                          tem_init_sparta, tem_destroy_sparta, &
    &                          tem_output_sparta
  use tem_logging_module,    only: logUnit, tem_logging_load_primary
  use tem_general_module,    only: tem_general_type, tem_start

  use aotus_module,          only: open_config_chunk
  use flu_binding,           only: flu_state

  implicit none

  ! MPI variables
  integer :: iError
  integer :: myrank
  integer :: nprocs
  integer :: comm
  integer :: stat(MPI_STATUS_SIZE)
  type(flu_state) :: conf
  type(tem_general_type) :: general
  character(len=120) :: outline

  logical :: OK = .false.
  logical :: correct = .false.

  ! main variables
  real(kind=rk), allocatable :: weight(:)
  integer :: myElems
  integer :: iRank
  integer(kind=long_k) :: offset
  type( tem_sparta_type ) :: sparta

  call tem_start('TREELM unit test', general)
  comm = general%proc%comm
  myrank = general%proc%rank
  nprocs = general%proc%comm_size
  if ( nprocs /= 5 ) stop

  ! Open the configuration file 
  call open_config_chunk(L = conf, chunk = trim(cubeconf))
  ! load and initialize logUnit
  call tem_logging_load_primary(conf = conf,  &
    &                           rank = myrank )

  allocate( weight(5) )
  myElems = 5
  select case ( myrank )
    case (0)
      weight = [ 5.0, 3.0, 1.0, 2.0, 1.0 ]
    case (1)
      weight = [ 4.0, 6.0, 1.0, 3.0, 2.0 ]
    case (2)
      weight = [ 1.0, 3.0, 1.0, 1.0, 1.0 ]
    case (3)
      weight = [ 1.0, 9.0, 1.0, 1.0, 1.0 ]
    case (4)
      weight = [ 1.0, 1.0, 1.0, 1.0, 1.0 ]
    case default
      stop
  end select

  call tem_init_sparta( sparta, nprocs )
  call tem_balance_sparta(weight, myrank, nprocs, comm, myElems, offset, &
    &                    sparta )
  call tem_output_sparta( sparta, logUnit(1) )

  write(outline,"(3(A,I2),A,F5.1)") "After balance, rank: ", myrank, &
    &                               ", myElems: ", myElems,          &
    &                               ", offset: ", offset,            &
    &                               ", my workload: ", sum(weight)

  if (myrank == 0 .and. myElems == 4 .and. offset == 0) then
    OK = .true.
  else if (myrank == 1 .and. myElems == 3 .and. offset == 4) then
    OK = .true.
  else if (myrank == 2 .and. myElems == 5 .and. offset == 7) then
    OK = .true.
  else if (myrank == 3 .and. myElems == 5 .and. offset == 12) then
    OK = .true.
  else if (myrank == 4 .and. myElems == 8 .and. offset == 17) then
    OK = .true.
  else
    OK = .false.
  end if

  if (myrank == 0) then

    write(*,*) trim(outline)

    do irank=1,nprocs-1
      call mpi_recv(outline, 120, MPI_CHARACTER, iRank, 42, comm, stat, iError)
      write(*,*) trim(outline)
    end do

  else

    call mpi_send(outline, 120, MPI_CHARACTER, 0, 42, comm, iError)

  end if

  call mpi_reduce( OK, correct, 1, mpi_logical, mpi_land, 0, comm, iError )

  call tem_destroy_sparta( sparta )

  call MPI_Barrier(comm, iError)

  if ( myrank == 0 ) then
    flush(6)
    if ( .not. correct ) then
      write(*,"(3(A,I0),A,F7.1)") "After  balance, rank: ", myrank, &
        &                      ", myElems: ", myElems, &
        &                      ", offset: ", offset, &
        &                      ", my workload: ", sum(weight)
      write(*,*) "FAILED"
    else
      write(*,*) "PASSED"
    end if
  end if

  deallocate( weight )

  call fin_env()

end program tem_sparta_test
