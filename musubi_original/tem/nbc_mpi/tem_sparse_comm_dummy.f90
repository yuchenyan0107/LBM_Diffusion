! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
!> A dummy module to avoid the need for non-blocking collectives.
!!
!! See the actual tem_sparse_comm_module for details, this module
!! merely serves as a placeholder to provide all the interfaces,
!! if MPI does not not support non-blocking collectives.
module tem_sparse_comm_module
  use aotus_module,          only: flu_state, aot_get_val
  use tem_aux_module,        only: tem_abort
  use tem_grow_array_module, only: grw_intArray_type, append, destroy
  use tem_logging_module,    only: logUnit

  implicit none

  private

  !> Data structure to store incoming requests.
  type serve_req_type
    !> Immediate buffer to receive data.
    integer :: recvdat

    !> Request handle for the opened receive.
    integer :: recv_req

    !> MPI communicator to do the communication in.
    integer :: comm

    !> Tag to use during MPI communications.
    integer :: tag

    !> An initial length to use for the growing arrays.
    integer :: inilen

    !> Source ranks we got requests from.
    type(grw_intArray_type) :: sources

    !> Growing array to store the sent data from source ranks.
    type(grw_intArray_type) :: src_dat
  end type serve_req_type

  public :: tem_sparse_alltoall_int
  public :: tem_sparse_comm_load

  !> Use the sparse all to all communication instead of the global MPI alltoall
  !!
  !! Defaults to false. It might be useful for large process counts,
  !! but requires nonblocking collectives from MPI-3.
  logical, save, public :: use_sparse_alltoall = .false.


contains


  !> Load the setting whether to use the sparse alltoall exhchange or not from
  !! the configuration.
  subroutine tem_sparse_comm_load(conf)
    !> Configuration handle to read the setting on the sparse alltoall from
    type(flu_state) :: conf

    integer :: iError

    ! Configure, if we are to write the Runtime info in tem_finalize
    call aot_get_val( L       = conf,                  &
      &               key     = 'use_sparse_alltoall', &
      &               val     = use_sparse_alltoall,   &
      &               ErrCode = iError,                &
      &               default = .false.                )

    if (use_sparse_alltoall) then
      write(logunit(1),*) 'ERROR: Requested sparse alltoall, but non-blocking'
      write(logunit(1),*) '       collectives were not available in the'
      write(logunit(1),*) '       MPI library, this executable was compiled'
      write(logunit(1),*) '       with!'
      write(logunit(1),*) ''
      write(logunit(1),*) 'If you want to make use of the sparse alltoall, you'
      write(logunit(1),*) 'need to use an MPI implementation supporting MPI-3.'
      write(logunit(1),*) ''
      write(logunit(1),*) 'Stopping now...'
      call tem_abort()
    end if

  end subroutine tem_sparse_comm_load


  !> Dummy for the sparse data exchange.
  !!
  subroutine tem_sparse_alltoall_int( targets, send_buffer, &
    &                                 sources, recv_buffer, &
    &                                 comm, tag             )

    !> List of target ranks to send an integer to.
    integer, intent(in) :: targets(:)

    !> Data to send to the respective target ranks. This array has to have the
    !! same ordering as targets.
    integer, intent(in) :: send_buffer(:)

    !> List of ranks we received data from (source ranks).
    !! The array will be allocated with a size according to the number of
    !! processes that send a request to this process.
    integer, intent(out), allocatable :: sources(:)

    !> Received data from the sources. The array has the same size and ordering
    !! as the sources array.
    integer, intent(out), allocatable :: recv_buffer(:)

    !> MPI Communicator to use for this data exchange.
    integer, intent(in) :: comm

    !> Tag to use in the communications. Defaults to 22.
    integer, intent(in), optional :: tag

    write(logunit(1),*) 'Sorry, no non-blocking collectives available.'
    write(logunit(1),*) 'Can not do a sparse_alltoall.'
    write(logunit(1),*) ''
    write(logunit(1),*) 'This should not have happened! Somebody forgot to call'
    write(logunit(1),*) 'tem_sparse_comm_load in the initialization.'
    write(logunit(1),*) 'Please make sure, your application actually calls this'
    write(logunit(1),*) 'config loading for the sparse communication.'
    write(logunit(1),*) ''
    write(logunit(1),*) 'Stopping now...'
    call tem_abort()

  end subroutine tem_sparse_alltoall_int


end module tem_sparse_comm_module
