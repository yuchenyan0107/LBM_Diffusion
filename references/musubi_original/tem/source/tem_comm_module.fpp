! Copyright (c) 2011-2013 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011-2013, 2016, 2018-2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2011-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012-2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2012, 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2012 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2012 Sathish Krishnan P S <s.krishnan@grs-sim.de>
! Copyright (c) 2013 Melven Zoellner <yameta@freenet.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2016-2017 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
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
! **************************************************************************** !!
!> This module provides the data structure for the communication
!! during the simulation.
!!
!! Several exchange methods are implemented. CoCo is heavily used here to allow
!! for a concise definition of exchanges for various data types.
!! The basic idea is to initialize the buffers, use them in exchanges and
!! finalize them when not needed anymore.
!! In the definition of the buffers, an array of positions in the original
!! linearized data array is used to describe the origin or target positions for
!! communicated data.
!!
!! @note If you introduce a new type to exchange, you will need to introduce
!!      appropriate CoCo copy statements everywhere.
!!
module tem_comm_module

  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer

  ! include treelm modules
  use mpi
  use env_module,            only: rk, rk_mpi, long_k, long_k_mpi
  use tem_aux_module,        only: tem_abort, check_mpi_error
  use tem_logging_module,    only: logUnit, tem_toStr
  use tem_grow_array_module, only: grw_intArray_type,     &
    &                              init, append, destroy
  use tem_dyn_array_module,  only: dyn_intArray_type, append, destroy, &
    &                              PositionOfVal
  use hvs_sizeof_module,     only: c_sizeof
  use mem_for_mpi_module,    only: alloc_mpif_mem, free_mpif_mem
  use tem_sparse_comm_module, only: use_sparse_alltoall, tem_sparse_alltoall_int

  ! include aotus modules
  use flu_binding,  only: flu_State
  use aotus_module, only: aot_get_val

  implicit none

  private

  public :: tem_communication_type
  public :: tem_commPattern_type
  public :: tem_comm_dumpType
  public :: tem_load_commPattern
  public :: tem_comm_init
  public :: tem_comm_count
  public :: tem_comm_createBuffer
  public :: tem_comm_alltoall_int
  public :: tem_comm_destroy

! Use CoCo here to make this buffer type generic for different types
?? text :: buffer_txt(tname, tstring)
  public :: tem_?tname?buffer_type

  !> Process-wise buffer for data of type ?tstring?
  !!
  !! This datatype is used to describe the exchange with a specific process, in
  !! case of explicit buffers it provides the memory for them.
  type tem_?tname?buffer_type

    !> Explicit buffer for data to be transferred
    ?tstring?, pointer :: val(:) => NULL()

    !> Explicit buffer in memory allocated by MPI
    type(c_ptr) :: mem_mpi

    !> position in the input vector from where to read the entries in
    !! val_?tname?
    !!
    !! @note JZ: in ATELES we use this to specify the positions of the cell
    !!          states that have to be sent.
    integer, allocatable :: pos(:)

    !> Number of values to exchange
    !!
    !! @note JZ: in ATELES this variable stores the number of coefficients we
    !!          transfer, i.e. number of cells to transfer times number of
    !!          degree of freedoms per cell times the number of scalar
    !!          variables.
    integer :: nVals

    !> Handle for the MPI-Datatype to describe the memory access,
    !! without explicit copying in the application.
    integer :: memIndexed
  end type tem_?tname?buffer_type
?? end text buffer_txt

! Let coco actually implement the buffer type given above for specific data
! types
?? copy :: buffer_txt(long, integer(kind=long_k))
?? copy :: buffer_txt(int, integer)
?? copy :: buffer_txt(real, real(kind=rk))

  ! ------------------------------------------------------------------------ !
  !> Description of communication data
  type tem_communication_type

    integer :: nProcs=0     !< amount of partitions to send to

    !> partition MPI rank
    integer,allocatable :: proc(:)

    !> How many data elements need to be exchanged with proc (per process).
    integer,allocatable :: nElemsProc(:)

    !> Request handle array
    integer,allocatable :: rqHandle(:)

    !> Data element positions in the actual arrays, used to built the pos
    !! information in the actual buffers (per process).
    type(grw_intArray_type), allocatable :: elemPos(:)

    !> declare communication buffers for each variable type
?? text :: buf_component_txt(tname)
    type( tem_?tname?buffer_type ), allocatable :: buf_?tname?(:)
?? end text buf_component_txt

?? copy :: buf_component_txt(long)
?? copy :: buf_component_txt(int)
?? copy :: buf_component_txt(real)
  end type tem_communication_type
  ! ------------------------------------------------------------------------ !


  ! ------------------------------------------------------------------------ !
  !> General description of the communication pattern to use.
  !!
  !! Depending on the chosen style, different exchange implementations are
  !! used. This data type provides the appropriate function pointers for
  !! initialization, finalization and exchange of the buffers.
  type tem_commPattern_type
    character(len=40) :: style

?? text :: exch_routine_txt(tname)
  procedure(tem_exchange_?tname?), nopass, pointer :: exchange_?tname?
  procedure(tem_commbuf_?tname?_init), nopass, pointer :: initBuf_?tname?
  procedure(tem_commbuf_?tname?_fin), nopass, pointer :: finBuf_?tname?
?? end text exch_routine_txt

?? copy :: exch_routine_txt(long)
?? copy :: exch_routine_txt(int)
?? copy :: exch_routine_txt(real)
  end type tem_commPattern_type
  ! ------------------------------------------------------------------------ !

! Definition of the abstract interfaces for the function pointers in the
! TEM_commPattern_type
?? text :: exchange_interface_txt(tname, tstring)
  abstract interface
    subroutine tem_exchange_?tname?( send, recv, state, message_flag, &
      &                              send_state, comm                 )
      import :: rk, long_k, tem_communication_type
      type(tem_communication_type), intent(inout) :: send, recv
      ?tstring?, intent(inout) :: state(*)
      integer, intent(in) :: message_flag
      ?tstring?, intent(in), optional :: send_state(*)
      !> MPI Communicator
      integer, intent(in) :: comm
    end subroutine tem_exchange_?tname?

    subroutine tem_commbuf_?tname?_init(me, pos, nVals)
      import :: tem_?tname?buffer_type
      type(tem_?tname?buffer_type), intent(inout) :: me
      integer, intent(in) :: nVals
      integer, intent(in) :: pos(nVals)
    end subroutine tem_commbuf_?tname?_init

    subroutine tem_commbuf_?tname?_fin(me)
      import :: tem_?tname?buffer_type
      type(tem_?tname?buffer_type), intent(inout) :: me
    end subroutine tem_commbuf_?tname?_fin
  end interface

?? end text exchange_interface_txt

?? copy :: exchange_interface_txt(long, integer(kind=long_k))
?? copy :: exchange_interface_txt(int, integer)
?? copy :: exchange_interface_txt(real, real(kind=rk))


contains


  ! ************************************************************************ !
  !> This subroutine loads the communication pattern from a Lua script
  !! and sets the exchange routine to be used accordingly.
  !!
  !! The variable read from the script is "commpattern".
  !! Several patterns are available:
  !! * isend_irecv - Use explicit buffers, copy first the outgoing ones, then
  !!      post irecvs and isends, wait on all, and copy the incoming buffers to
  !!      their final locations (default).
  !! * isend_irecv_mpimem - Same as isend_irecv, but use memory that is
  !!      allocatd by MPI_Alloc_mem for the buffers.
  !! * isend_irecv_overlap - Similar to isend_irecv, but directly post sends,
  !!      after filling the outgoing buffer to each process, wait on any
  !!      incoming messages and only wait on sends, after everything is copied
  !!      into the final location.
  !! * overlap_mpimem - Same as isend_irecv_overlap, but use memory that is
  !!      allocated by MPI_Alloc_mem for the buffers.
  !! * typed_isend_irecv - Instead of copying the memory around, define a
  !!      indexed MPI datatype and use that in the exchange.
  !! * gathered_type - Similar to typed_isend_irecv, but with minimal number of
  !!      blocks in the indexed type, by tracking only contiguous blocks in the
  !!      memory layout.
  !!
  !! Instead of reading the style from a configuration script, it can also be
  !! directly set by the caller. If nothing is specified by the caller, the
  !! style will default to isend_irecv.
  !!
  !! Usage:
  !! ~~~~~~~~~~~~~~~~~~~~~{.lua}
  !! commpattern = 'isend_irecv'
  !! ~~~~~~~~~~~~~~~~~~~~~
  !!
  subroutine tem_load_commPattern( me, conf, style )
    ! -------------------------------------------------------------------- !
    !> commpattern to set
    type(tem_commPattern_type), intent(out) :: me
    !> handle to the Lua script
    type(flu_State), optional :: conf
    !> optional communication style
    character(len=*), intent(in), optional :: style
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    write(logUnit(1),*)"Loading communication pattern:"

    if (present(conf)) then
      ! If a configuration is given, this trumps any other setting.
      ! Defaults to isend_irecv.
      call aot_get_val( L       = conf,          &
        &               key     = 'commpattern', &
        &               val     = me%style,      &
        &               ErrCode = iError,        &
        &               default = 'isend_irecv'  )
    else if (present(style)) then
      ! If a style is given directly by the caller, use that one.
      me%style = style
    else
      ! Default to isend_irecv if nothing provided by the caller.
      me%style = 'isend_irecv'
    end if

?? text :: exch_pointer_txt(routinename, initname, finname, tname)
      me%exchange_?tname? => ?routinename?_?tname?
      me%initBuf_?tname? => tem_commbuf_?tname?_?initname?
      me%finBuf_?tname? => tem_commbuf_?tname?_?finname?
?? end text exch_pointer_txt

    select case(trim(me%style))
    case ('isend_irecv')
?? copy :: exch_pointer_txt(comm_isend_irecv, fillPos, finPos, long)
?? copy :: exch_pointer_txt(comm_isend_irecv, fillPos, finPos, int)
?? copy :: exch_pointer_txt(comm_isend_irecv, fillPos, finPos, real)

    case ('isend_irecv_overlap')
?? copy :: exch_pointer_txt(comm_isend_irecv_overlap, fillPos, finPos, long)
?? copy :: exch_pointer_txt(comm_isend_irecv_overlap, fillPos, finPos, int)
?? copy :: exch_pointer_txt(comm_isend_irecv_overlap, fillPos, finPos, real)

    case ('typed_isend_irecv')
?? copy :: exch_pointer_txt(comm_typed_isend_irecv, fillIndexed, finTyped, long)
?? copy :: exch_pointer_txt(comm_typed_isend_irecv, fillIndexed, finTyped, int)
?? copy :: exch_pointer_txt(comm_typed_isend_irecv, fillIndexed, finTyped, real)

    case ('gathered_type')
?? copy :: exch_pointer_txt(comm_typed_isend_irecv, gatherIndexed, finTyped, long)
?? copy :: exch_pointer_txt(comm_typed_isend_irecv, gatherIndexed, finTyped, int)
?? copy :: exch_pointer_txt(comm_typed_isend_irecv, gatherIndexed, finTyped, real)

    case ('isend_irecv_mpimem')
?? copy :: exch_pointer_txt(comm_isend_irecv, fillmpimem, finmpimem, long)
?? copy :: exch_pointer_txt(comm_isend_irecv, fillmpimem, finmpimem, int)
?? copy :: exch_pointer_txt(comm_isend_irecv, fillmpimem, finmpimem, real)

    case ('overlap_mpimem')
?? copy :: exch_pointer_txt(comm_isend_irecv_overlap, fillmpimem, finmpimem, long)
?? copy :: exch_pointer_txt(comm_isend_irecv_overlap, fillmpimem, finmpimem, int)
?? copy :: exch_pointer_txt(comm_isend_irecv_overlap, fillmpimem, finmpimem, real)

    case default
      write(logUnit(1),*) "ERROR, unknown commpattern: "//trim(me%style)
      write(logUnit(1),*) "available are: "
      write(logUnit(1),*) "* isend_irecv"
      write(logUnit(1),*) "* isend_irecv_overlap"
      write(logUnit(1),*) "* typed_isend_irecv"
      write(logUnit(1),*) "* gathered_type"
      write(logUnit(1),*) "* isend_irecv_mpimem"
      write(logUnit(1),*) "* overlap_mpimem"
      call tem_abort()

    end select

    write(logUnit(1),*) trim(me%style)

  end subroutine tem_load_commPattern
  ! ************************************************************************ !

! The following template defines the routines for various communication
! patterns, if you want to add a new pattern, you have to define an init,
! finalize and exchange routine in this template.
! To add a new data type, just add the appropriate CoCo copy line after this
! template.
?? text :: exchange_impl_txt(tname, tstring, vartype_mpi )
  ! ************************************************************************ !
  !> Fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer.
  !!
  subroutine tem_commbuf_?tname?_fillPos( me, pos, nVals )
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    integer, intent(in) :: nVals
    integer, intent(in) :: pos(nVals)
    ! -------------------------------------------------------------------- !

    me%nVals = nVals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nVals))
    me%pos = pos

    if ( associated(me%val) ) deallocate(me%val)
    allocate(me%val(nVals))

  end subroutine tem_commbuf_?tname?_fillPos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Free the communication buffer allocated by the fillPos routine again.
  !!
  subroutine tem_commbuf_?tname?_finPos(me)
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nVals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) deallocate(me%val)

  end subroutine tem_commbuf_?tname?_finPos
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Fill the positions that describe how the data in the state
  !! vector relates to the entries in the buffer use memory that
  !! is allocated by MPI for the buffer.
  !!
  subroutine tem_commbuf_?tname?_fillmpimem( me, pos, nVals )
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    integer, intent(in) :: nVals
    integer, intent(in) :: pos(nVals)
    ! -------------------------------------------------------------------- !
    ?tstring? :: typesample
    integer(kind=MPI_ADDRESS_KIND) :: typelen
    ! -------------------------------------------------------------------- !

    typelen = int(c_sizeof(typesample), kind=MPI_ADDRESS_KIND)

    me%nVals = nVals

    if ( allocated(me%pos) ) deallocate(me%pos)
    allocate(me%pos(nVals))
    me%pos = pos

    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if
    call alloc_mpif_mem( asize   = nVals*typelen, &
      &                  baseptr = me%mem_mpi     )
    call c_f_pointer(me%mem_mpi, me%val, [me%nVals])

  end subroutine tem_commbuf_?tname?_fillmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Free the communication buffer allocated by the fillmpimem routine again.
  !!
  subroutine tem_commbuf_?tname?_finmpimem(me)
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !

    me%nVals = 0
    if ( allocated(me%pos) ) deallocate(me%pos)
    if ( associated(me%val) ) then
      nullify(me%val)
      call free_mpif_mem(me%mem_mpi)
    end if

  end subroutine tem_commbuf_?tname?_finmpimem
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Fill the indexed MPI datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  subroutine tem_commbuf_?tname?_fillIndexed(me, pos, nVals)
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    integer, intent(in) :: nVals
    integer, intent(in) :: pos(nVals)
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    me%nVals = nVals
    ! CALL MPI_TYPE_CREATE_INDEXED_BLOCK(COUNT, BLOCKLENGTH, &
    !   &                   ARRAY_OF_DISPLACEMENTS, &
    !   &                   OLDTYPE, NEWTYPE, IERROR)
    call MPI_Type_create_indexed_block( nVals, 1, pos - 1,                   &
      &                                 ?vartype_mpi?, me%memIndexed, iError )
    call check_mpi_error(iError,'create indexed block in tem_commbuf_?tname?_fillIndexed')
    call MPI_Type_commit(me%memIndexed, iError)
    call check_mpi_error(iError,'commit memIndexed in tem_commbuf_?tname?_fillIndexed')

  end subroutine tem_commbuf_?tname?_fillIndexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Free the communication buffer allocated by the fillIndexed routine again.
  !!
  subroutine tem_commbuf_?tname?_finTyped(me)
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    ! -------------------------------------------------------------------- !
    integer :: iError
    ! -------------------------------------------------------------------- !

    me%nVals = 0
    call MPI_Type_free(me%memIndexed, iError)
    call check_mpi_error(iError,'free memIndexed in tem_commbuf_?tname?_finTyped')

  end subroutine tem_commbuf_?tname?_finTyped
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Gather the indexed MPI datatype, which describes how the data in the state
  !! vector relates to the entries in the buffer.
  !! In contrast to the simple indexed type above, we try to minimize the number
  !! of blocks here, and gather contiguous blocks of memory together.
  !!
  subroutine tem_commbuf_?tname?_gatherIndexed( me, pos, nVals )
    ! -------------------------------------------------------------------- !
    type(tem_?tname?buffer_type), intent(inout) :: me
    integer, intent(in) :: nVals
    integer, intent(in) :: pos(nVals)
    ! -------------------------------------------------------------------- !
    type(grw_IntArray_type) :: blocklength
    type(grw_IntArray_type) :: displ
    integer :: iVal, counter
    integer :: iError
    ! -------------------------------------------------------------------- !

    me%nVals = nVals

    ! Initialize growing arrays, a kB should be fine to start with...
    call init(blocklength, 256)
    call init(displ, 256)

    if (nVals > 0) then

      ! Start with the displacement of the first entry in the list
      call append(displ, pos(1)-1)
      counter = 1

      do iVal=2,nVals
        if (pos(iVal) == pos(iVal-1)+1) then
          ! Contiguous memory location following the previous one, increase the
          ! the blocklength.
          counter = counter + 1
        else
          ! New block encountered, record the block found so far
          call append(blocklength, counter)

          ! Start new block
          call append(displ, pos(iVal)-1)
          counter = 1
        end if
      end do

      ! Finish the last block, by recording its found length:
      call append(blocklength, counter)

    end if

    ! CALL MPI_TYPE_INDEXED(COUNT, ARRAY_OF_BLOCKLENGTHS, &
    !   &                   ARRAY_OF_DISPLACEMENTS, OLDTYPE, NEWTYPE, IERROR)
    call MPI_Type_indexed( displ%nVals, blocklength%Val, displ%Val, &
      &                    ?vartype_mpi?, me%memIndexed, iError     )
    call check_mpi_error(iError,'type indexed in tem_commbuf_?tname?_gatherIndexed')
    call MPI_Type_commit(me%memIndexed, iError)
    call check_mpi_error(iError,'commit memIndexed in tem_commbuf_?tname?_gatherIndexed')

    call destroy(displ)
    call destroy(blocklength)

  end subroutine tem_commbuf_?tname?_gatherIndexed
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Exchange the communication buffers
  !! with a non-blocking MPI communication using preposted iRecv and iSend
  !! with a waitall
  !!
  subroutine comm_isend_irecv_?tname?( send, recv, state, message_flag, &
    &                                  send_state, comm                 )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send
    type( tem_communication_type ), intent(inout) :: recv
    ?tstring?, intent(inout) :: state(*)  !< State vector to update
    integer, intent(in) :: message_flag
    ?tstring?, intent(in), optional :: send_state(*)  !< Data to send
    !> MPI Communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! request handle for messages
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nProcs + send%nProcs )
    integer :: status( mpi_status_size, max(recv%nProcs, send%nProcs) )
    integer :: iErr             !< error flag
    integer :: iProc, iVal
    integer :: nSendVals, nRecvVals
    ! -------------------------------------------------------------------- !

    if (present(send_state)) then
      do iProc = 1, send%nProcs
        ! Fill communication message
        nSendVals = send%buf_?tname?( iProc )%nVals
        !$NEC ivdep
        do iVal = 1, nSendVals
          send%buf_?tname?( iProc )%val( iVal )                    &
            &  = send_state( send%buf_?tname?( iProc )%pos( iVal ) )
        end do
      end do
    else
      do iProc = 1, send%nProcs
        ! Fill communication message
        nSendVals = send%buf_?tname?( iProc )%nVals
        !$NEC ivdep
        do iVal = 1, nSendVals
          send%buf_?tname?( iProc )%val( iVal )               &
            &  = state( send%buf_?tname?( iProc )%pos( iVal ) )
        end do
      end do
    end if

    do iProc = 1, recv%nProcs
      ! Start receive communications
      call mpi_irecv(                           &
       &      recv%buf_?tname?( iProc )%val,    & ! me
       &      recv%buf_?tname?( iProc )%nVals,  & ! me size
       &      ?vartype_mpi?,                    & ! data type
       &      recv%proc(iProc),                 & ! target me
       &      message_flag,                     & ! flag
       &      comm,                             & ! communicator
       &      recv%rqHandle(iProc),             & ! handle
       &      iErr )                              ! error status
    enddo

    !  Start the sending communications
    do iProc = 1, send%nProcs
      call mpi_isend(                           &
       &      send%buf_?tname?( iProc )%val,    & ! buffer
       &      send%buf_?tname?( iProc )%nVals,  & ! count
       &      ?vartype_mpi?,                    & ! data type
       &      send%proc(iProc),                 & ! target
       &      message_flag,                     & ! tag
       &      comm,                             & ! communicator
       &      send%rqHandle( iProc ),           & ! handle
       &      iErr )                              ! error status
    enddo ! iProc

    ! Wait for receive buffer to be ready
    if ( recv%nProcs /= 0 ) then
      call mpi_waitall(recv%nProcs,               & ! count
        &              recv%rqHandle,             & ! request handles
        &              status,                    & ! mpi status
        &              iErr )                       ! error status
    end if

    ! Now values from recv me can be copied to the actual state array
    do iProc = 1, recv%nProcs
      nRecvVals = recv%buf_?tname?( iProc )%nVals
      !$NEC ivdep
      do iVal = 1, nRecvVals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_?tname?( iProc )%pos( iVal ) ) &
          &  = recv%buf_?tname?( iProc )%val( iVal )
      end do
    end do

    ! Wait for send buffer to be ready
    if ( send%nProcs /= 0 ) then
      call mpi_waitall(send%nProcs,   & ! count
        &              send%rqHandle, & ! request handles
        &              status,        & ! mpi status
        &              iErr )           ! error status
    end if

  end subroutine comm_isend_irecv_?tname?
  ! ************************************************************************ !


  ! ************************************************************************ !
  subroutine comm_isend_irecv_overlap_?tname?( send, recv, state,             &
    &                                          message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type), intent(inout) :: send, recv
    ?tstring?, intent(inout) :: state(*)
    integer, intent(in) :: message_flag
    ?tstring?, intent(in), optional :: send_state(*)  !< Data to send
    !> MPI Communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    integer :: iProc, iVal  ! counter for neigbor mees
    integer :: finProc
    ! integer :: rq_handle( recv%nProcs + send%nProcs )
    integer :: status( mpi_status_size, max(send%nProcs,1) )
    integer :: iErr ! error flag
    integer :: nSendVals, nRecvVals
    ! -------------------------------------------------------------------- !

    ! Pre-post iRecv
    do iProc = 1, recv%nProcs
      call mpi_irecv(                          &
       &      recv%buf_?tname?( iProc )%val,   & ! buffer
       &      recv%buf_?tname?( iProc )%nVals, & ! count
       &      ?vartype_mpi?,                   & ! data type
       &      recv%proc(iProc),                & ! target me
       &      message_flag,                    & ! tag
       &      comm,                            & ! communicator
       &      recv%rqHandle(iProc),            & ! handle
       &      iErr                             ) ! error status
    end do

    !> Fill send buffers and start Sending
    do iProc = 1, send%nProcs
      nSendVals = send%buf_?tname?( iProc )%nVals
      if (present(send_state)) then
        !$NEC ivdep
        do iVal = 1, nSendVals
          send%buf_?tname?( iProc )%val( iVal )                    &
            &  = send_state( send%buf_?tname?( iProc )%pos( iVal ) )
        end do
      else
        !$NEC ivdep
        do iVal = 1, nSendVals
          send%buf_?tname?( iProc )%val( iVal )               &
            &  = state( send%buf_?tname?( iProc )%pos( iVal ) )
        end do
      end if
      call mpi_isend(                               &
       &      send%buf_?tname?( iProc )%val,        & ! buffer
       &      send%buf_?tname?( iProc )%nVals,      & ! count
       &      ?vartype_mpi?,                        & ! data type
       &      send%proc(iProc),                     & ! target
       &      message_flag,                         & ! tag
       &      comm,                                 & ! comm
       &      send%rqHandle( iProc ),               & ! handle
       &      iErr )                                  ! error status
    end do

    do iProc = 1, recv%nProcs
      ! Wait for any of receive buffer to be ready
      call mpi_waitany(   &
        &    recv%nProcs, & ! count
        &    recv%rqHandle,   & ! request handles
        &    finProc,     & ! process that finished
        &    status(:,1), & ! mpi status
        &    iErr )         ! error status
      nRecvVals = recv%buf_?tname?( finProc )%nVals
      !$NEC ivdep
      do iVal = 1, nRecvVals
        ! write the values from the recv me to the state array
        ! to positions given in recvme%pos
        state( recv%buf_?tname?( finProc )%pos( iVal ) ) &
          &  = recv%buf_?tname?( finProc )%val( iVal )
      end do
    end do

    if (send%nProcs > 0) then
      ! Wait for send buffer to be ready
      call mpi_waitall(                  &
        &    send%nProcs,                & ! total number of comm.'s to wait for
        &    send%rqHandle,              & ! request handles
        &    status,                     & ! mpi status
        &    iErr )                        ! error status
    end if

  end subroutine comm_isend_irecv_overlap_?tname?
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Exchange the communication mes
  !! with a non-blocking MPI communication using preposted iRecv and iSend
  !! with a waitall
  !!
  subroutine comm_typed_isend_irecv_?tname?( send, recv, state,             &
    &                                        message_flag, send_state, comm )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: send, recv
    ?tstring?, intent(inout) :: state(*)  !< Current state vector
    integer, intent(in) :: message_flag
    ?tstring?, intent(in), optional :: send_state(*)  !< Data to send
    !> MPI Communicator
    integer, intent(in) :: comm
    ! -------------------------------------------------------------------- !
    ! @todo request handle array could exist during complete code runtime
    ! integer :: rq_handle( recv%nProcs + send%nProcs )
    integer :: status( mpi_status_size, max(recv%nProcs, send%nProcs) )
    integer :: iErr             ! error flag
    integer :: iProc
    ! -------------------------------------------------------------------- !

    !> Values for send me must have been copied from the actual state array
    do iProc = 1, recv%nProcs
      !> Start receive communications
      call mpi_irecv(                                     &
       &              state,                              & ! buffer
       &              1,                                  & ! count
       &              recv%buf_?tname?(iProc)%memIndexed, & ! type
       &              recv%proc(iProc),                   & ! source
       &              message_flag,                       & ! tag
       &              comm,                               & ! comm
       &              recv%rqHandle(iProc),               & ! request handle
       &              iErr )

    end do

    !>  Start the sending communications
    if (present(send_state)) then
      do iProc = 1, send%nProcs
        call mpi_isend(                                     &
          &             send_state,                         & ! buffer
          &             1,                                  & ! count
          &             send%buf_?tname?(iProc)%memIndexed, & ! type
          &             send%proc(iProc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqHandle( iProc ),             & ! handle
          &             iErr )
      end do !< iProc
    else
      do iProc = 1, send%nProcs
        call mpi_isend(                                     &
          &             state,                              & ! buffer
          &             1,                                  & ! count
          &             send%buf_?tname?(iProc)%memIndexed, & ! type
          &             send%proc(iProc),                   & ! target
          &             message_flag,                       & ! tag
          &             comm,                               & ! comm
          &             send%rqHandle( iProc ),             & ! handle
          &             iErr )
      end do !< iProc
    end if

    !> Wait for above communications to complete
    if ( recv%nProcs /= 0 ) then
      call mpi_waitall( recv%nProcs,     & ! count
        &               recv%rqHandle,   & ! request handles
        &               status,          & ! statuses
        &               iErr )
    end if
    if ( send%nProcs /= 0 ) then
      call mpi_waitall( send%nProcs,     & ! count
        &               send%rqHandle,   & ! request handles
        &               status,          & ! statuses
        &               iErr )
    end if

  end subroutine comm_typed_isend_irecv_?tname?
  ! ************************************************************************ !
?? end text exchange_impl_txt

?? copy :: exchange_impl_txt(long, integer(kind=long_k), long_k_mpi)
?? copy :: exchange_impl_txt(int, integer, mpi_integer)
?? copy :: exchange_impl_txt(real, real(kind=rk), rk_mpi )


  ! ************************************************************************ !
  !> Write communication type data to nUnit (debugging routine)
  !!
  subroutine tem_comm_dumpType(me, nUnit)
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(in) :: me
    integer, intent(in) :: nUnit
    ! -------------------------------------------------------------------- !

    write(nUnit,*)        '----   Communication type    -----------------------'
    write(nUnit,"(A,I0)") '       nProcs: ', me%nProcs
    write(nUnit,"(A   )") '         proc: '//trim(tem_toStr( me%proc, ',' ))
    write(nUnit,"(A   )") '   nElemsProc: ' &
      &                   // trim(tem_toStr( me%nElemsProc, ',' ))
    write(nUnit,*)        '---------------------------------------------'

  end subroutine tem_comm_dumpType
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Allocate tem_communication_type and its variables
  !!
  subroutine tem_comm_init( me, nProcs )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: me
    integer, intent(in) :: nProcs
    ! -------------------------------------------------------------------- !

    if ( allocated(me%proc) )       deallocate( me%proc )
    if ( allocated(me%nElemsProc) ) deallocate( me%nElemsProc )
    if ( allocated(me%elemPos) )    deallocate( me%elemPos )
    if ( allocated(me%rqHandle) )   deallocate( me%rqHandle )
    if ( allocated(me%buf_long) )   deallocate( me%buf_long )
    if ( allocated(me%buf_int) )    deallocate( me%buf_int )
    if ( allocated(me%buf_real) )   deallocate( me%buf_real )

    me%nProcs = nProcs
    allocate( me%proc      ( nProcs ) )
    allocate( me%nElemsProc( nProcs ) )
    allocate( me%elemPos   ( nProcs ) )
    allocate( me%rqHandle  ( nProcs ) )

    allocate( me%buf_long  ( nProcs ) )
    allocate( me%buf_int   ( nProcs ) )
    allocate( me%buf_real  ( nProcs ) )

  end subroutine tem_comm_init
  ! ************************************************************************ !


  ! ************************************************************************ !
  !> Allocate tem_communication_type and its variables
  !!
  subroutine tem_comm_count( me, comm_size, nHalos )
    ! -------------------------------------------------------------------- !
    type( tem_communication_type ), intent(inout) :: me
    !> communicator size
    integer, intent(in) :: comm_size
    !> number of halos for each other processes
    integer, intent(in) :: nHalos( comm_size )
    ! -------------------------------------------------------------------- !
    integer :: iPartner, iProc
    ! -------------------------------------------------------------------- !

    iPartner = 0
    do iProc = 1, comm_size
      if( nHalos( iProc ) > 0 ) then
        ! Store the processes numbers to receive from
        iPartner = iPartner + 1
        me%nElemsProc( iPartner ) = nHalos( iProc )
        me%proc( iPartner )       = iProc - 1
      end if
    end do

  end subroutine tem_comm_count
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> Routine to build communication buffer using elemRanks.
  !! This routine can be used only if all elements need to be communicated
  !! but they need process-wise seperation.
  !! Uses nScalars to get position in the value array to communicate.
  !! For send buffer: elemRanks contains target ranks to send data to
  !! For recv buffer: elemRanks contains source ranks to recv data from
  subroutine tem_comm_createBuffer( commBuffer, nScalars, nElems, elemRanks )
    ! -------------------------------------------------------------------------!
    !> send or recv communication buffer to be created
    type(tem_communication_type), intent(out) :: commBuffer
    !> Number of scalars per element
    integer, intent(in) :: nScalars
    !> Total number of elements or points to communicate
    integer, intent(in) :: nElems
    !> Target or source rank for each element or point
    integer, intent(in) :: elemRanks(nElems)
    ! -------------------------------------------------------------------------!
    integer :: iElem, iProc, iVar, counter, pntPos
    type(dyn_intArray_type) :: partnerProc
    integer, allocatable :: pos(:)
    ! -------------------------------------------------------------------------!
    ! Create dynamic array of rank ids to communicate to initialize commBuffer
    do iElem = 1, nElems
      call append( me  = partnerProc,    &
        &          val = elemRanks(iElem) )
    end do

    ! Initialize commBuffer
    call tem_comm_init( me     = commBuffer,       &
      &                 nProcs = partnerProc%nVals )
    ! Store rank id to communicate
    do iProc = 1, partnerProc%nVals
      commBuffer%proc(iProc) = partnerProc%val(iProc)
    end do

    ! Create map from commBuffer to elem array
    do iElem = 1, nElems
      iProc = PositionOfVal( me  = partnerProc,     &
        &                    val = elemRanks(iElem) )
      call append( me  = commBuffer%elemPos(iProc), &
        &          val = iElem                       )
    end do
    call destroy(partnerProc)

    commBuffer%nElemsProc(:) = commBuffer%elemPos(:)%nVals

    ! Get position in state array to send or recv data
    allocate( pos( nScalars * maxVal( commBuffer%nElemsProc(:) ) ) )

    ! Assign comm buffer positions
    do iProc = 1, commBuffer%nProcs
      counter = 0

      ! loop of nElems per proc
      do iElem = 1, commBuffer%nElemsProc(iProc)
        ! position of this proc point in the point array
        pntPos = commBuffer%elemPos(iProc)%val(iElem)
        do iVar = 1, nScalars
          counter = counter + 1
          ! position in evalVal array which has size: nElems*nScalars
          pos(counter) = (pntPos-1)*nScalars + iVar
        end do !iVar
      end do !iElem
      ! copy position array to me%pos, allocate me%val array
      call tem_commbuf_real_fillPos( me    = commBuffer%buf_real(iProc), &
        &                            pos   = pos,                        &
        &                            nVals = counter                     )
    end do !iProc
    deallocate(pos)

  end subroutine tem_comm_createBuffer
  ! ************************************************************************ !

  ! ************************************************************************ !
  !> All to all exchange of a single integer.
  !!
  !! This is a wrapper around the sparse alltoall implementation and overcome
  !! the lack of non-blocking collectives on some systems.
  subroutine tem_comm_alltoall_int( targets, send_buffer, &
    &                               sources, recv_buffer, &
    &                               comm, tag             )

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

    integer :: nProcs
    integer :: nSources
    integer :: iProc, iSource
    integer :: iError
    integer, allocatable :: buf(:)

    if (use_sparse_alltoall) then

      call tem_sparse_alltoall_int( targets     = targets,     &
        &                           send_buffer = send_buffer, &
        &                           sources     = sources,     &
        &                           recv_buffer = recv_buffer, &
        &                           comm        = comm,        &
        &                           tag         = tag          )

    else

      call MPI_Comm_size(comm, nProcs, iError)
      allocate(buf(0:nProcs-1))
      buf = 0
      buf(targets(:)) = send_buffer
      call MPI_Alltoall( MPI_IN_PLACE, 1, MPI_INTEGER,     &
        &                buf, 1, MPI_INTEGER, comm, iError )
      nSources = count(buf/=0)
      allocate(sources(nSources))
      allocate(recv_buffer(nSources))
      recv_buffer = 0
      iSource = 1
      do iProc=0,nProcs-1
        if (buf(iProc) /= 0) then
          sources(iSource) = iProc
          recv_buffer(iSource) = buf(iProc)
          iSource = iSource + 1
        end if
      end do
      deallocate(buf)
    end if

  end subroutine tem_comm_alltoall_int
  ! *************************************************************************** !

  ! *************************************************************************** !
  subroutine tem_comm_destroy( me, commPattern )
    ! -------------------------------------------------------------------------!
    !> communication type to be destroyed
    type(tem_communication_type), intent(inout) :: me
    !> Communication pattern
    type(tem_commPattern_type), intent(in) :: commPattern
    ! -------------------------------------------------------------------------!
    integer :: iProc
    ! -------------------------------------------------------------------------!

    do iProc = 1, me%nProcs
      call commPattern%finbuf_real( me%buf_real(iProc) )
      call commPattern%finbuf_long( me%buf_long(iProc) )
      call commPattern%finbuf_int(  me%buf_int(iProc)  )
    end do

    deallocate( me%buf_real )
    deallocate( me%buf_long )
    deallocate( me%buf_int  )

    deallocate( me%proc )
    deallocate( me%nElemsProc )
    do iProc = 1, me%nProcs
      call destroy( me%elemPos(iProc) )
    end do
    deallocate( me%elemPos )
    deallocate( me%rqHandle )

  end subroutine tem_comm_destroy
  ! *************************************************************************** !

end module tem_comm_module
! **************************************************************************** !
