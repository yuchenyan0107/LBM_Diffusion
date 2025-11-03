! Copyright (c) 2016-2017 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
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
!> A module to provide sparse communication.
!!
!! Sometimes there is information that is required by other processes, but only
!! the requesting process knows about this.
!! As the requested processes are unaware of which ranks might talk to them,
!! they potentially have to listen for every other rank.
!!
!! A traditional resolution to this would be offered by collective
!! communications. However, this quickly gets too expensive if there are only
!! few actual communication partners but in total many computing ranks.
!!
!! For this scenario, this module provides a solution with point-to-point
!! communications and a non-blocking barrier to overcome the termination
!! problem.
module tem_sparse_comm_module
  use mpi
  use aotus_module, only: flu_state, aot_get_val
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
      write(logunit(2),*) 'Using SPARSE alltoall algorithm (requires MPI3)'
    else
      write(logunit(2),*) 'Using global alltoall'
    end if

  end subroutine tem_sparse_comm_load


  !> A sparse data exchange, where potentially each process may talk to any
  !! other process.
  !!
  !! This is mainly useful in the case where there are only few actual
  !! communication partners in a large communicator.
  !! Each process provides the ranks it wants to talk with and the data to send
  !! to them.
  !! The data to be received is left open and automatically determined by this
  !! exchange.
  !!
  !! This routine exchanges a single integer between processes, which might be
  !! used to initiate further communication later on in a traditional point to
  !! point approach.
  !!
  !! This routine is meant to offer a drop-in replacement of the traditional
  !! global MPI all-to-all routine.
  !!
  !! Each rank may only appear once in the targets list.
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

    integer :: nTargets
    integer :: iTarget
    type(serve_req_type) :: me
    integer :: send_req(size(targets))
    integer :: send_stat(MPI_STATUS_SIZE, size(targets))
    integer :: bar_req
    integer :: bar_stat(MPI_STATUS_SIZE)
    integer :: iError
    logical :: allSent, allProcsDone

    nTargets = size(targets)

    call init_serving( me     = me, &
      &                comm   = comm,           &
      &                tag    = tag,            &
      &                inilen = max(nTargets,1) )

    ! Send data to all targets
    do iTarget=1,nTargets
      ! MPI_ISSEND(BUF, COUNT, DATATYPE, DEST, TAG, COMM, REQUEST, IERROR)
      call MPI_iSsend(send_buffer(iTarget), 1, MPI_INTEGER, targets(iTarget), &
        &             me%tag, comm, send_req(iTarget), iError     )
      call serve_request(me)
    end do

    ! Wait on sending to complete.
    do
      ! MPI_TESTALL(COUNT, ARRAY_OF_REQUESTS, FLAG, ARRAY_OF_STATUSES, IERROR)
      call MPI_Testall(nTargets, send_req, allSent, send_stat, iError)
      call serve_request(me)
      if (allSent) EXIT
    end do

    ! Signal termination of the algorithm for this process.
    ! MPI_IBARRIER(COMM, REQUEST, IERROR)
    call MPI_iBarrier(comm, bar_req, iError)

    ! Continue serving requests until all processes signalled the completion
    ! of their sending part.
    do
      call serve_request(me)

      ! MPI_TEST(REQUEST, FLAG, STATUS, IERROR)
      call MPI_Test(bar_req, allProcsDone, bar_stat, iError)
      if (allProcsDone) EXIT
    end do

    ! Exchange complete, finalize and return the source ranks with their
    ! sent data to the caller.
    allocate(sources(me%sources%nVals))
    allocate(recv_buffer(me%src_dat%nVals))
    recv_buffer=0
    if (me%sources%nVals > 0) then
      sources     = me%sources%val(:me%sources%nVals)
      recv_buffer = me%src_dat%val(:me%src_dat%nVals)
    end if

    call destroy(me%sources)
    call destroy(me%src_dat)

    ! Ensure that all processes had the opportunity to shut down their pending
    ! receives to avoid any confusion with subsequent calls of this sparse
    ! alltoall exchange.
    call MPI_Barrier(comm, iError)

  end subroutine tem_sparse_alltoall_int


  !> Initialising the serving part.
  !!
  !! We need to listen for messages from any source rank with the tag 22.
  subroutine init_serving(me, comm, inilen, tag)
    !> Serving data structure to initialize.
    type(serve_req_type), intent(out) :: me

    !> MPI communicator to use for the data exchange.
    integer, intent(in) :: comm

    !> An initial length for the list of sources.
    integer, intent(in) :: inilen

    !> Tag to use during MPI communications
    integer, intent(in), optional :: tag

    me%comm = comm
    me%inilen = inilen

    if (present(tag)) then
      me%tag = tag
    else
      me%tag = 22
    end if

  end subroutine init_serving


  !> Serve any incoming requests.
  !!
  !! Check for incoming messages and store them if there has been one received.
  !! This needs to be repeated until we can be sure that all processes finished
  !! sending data.
  subroutine serve_request(me)
    !> Serving data structure to cater now.
    type(serve_req_type), intent(inout) :: me

    integer :: iError
    integer :: source
    integer :: recv_stat(MPI_STATUS_SIZE)
    logical :: got_request

    ! Check possible incoming message
    ! MPI_IPROBE(SOURCE, TAG, COMM, FLAG, STATUS, IERROR)
    call MPI_Iprobe(MPI_ANY_SOURCE, me%tag, me%comm, got_request, recv_stat, iError)

    if (got_request) then
      ! Store received request.
      source = recv_stat(MPI_SOURCE)

      ! Receive meesage
      ! MPI_RECV(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
      call MPI_Recv( me%recvdat, 1, MPI_INTEGER, source, &
        &            me%tag, me%comm, recv_stat, iError  )

      call append( me     = me%sources, &
        &          val    = source,     &
        &          length = me%inilen   )
      call append( me     = me%src_dat, &
        &          val    = me%recvdat, &
        &          length = me%inilen   )

    end if

  end subroutine serve_request

end module tem_sparse_comm_module
