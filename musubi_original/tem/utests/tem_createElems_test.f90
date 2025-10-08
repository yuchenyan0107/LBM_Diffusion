! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2015 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
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
! This program tries to evaluate the memory demand per element
!
program elem_memAlloc

  ! include treelm modules
  use env_module,         only: long_k
  use tem_debug_module,   only: tem_reportStatus, tem_debug_type
  use tem_stencil_module, only: tem_stencilElement_type
  use tem_element_module, only: tem_element_type, getSize, getReqSize,init,    &
    &                           append, eT_fluid, destroy
  use tem_logging_module, only: tem_logging_init

  implicit none
  ! ----------------------------------------------------------------------------
  type(tem_element_type) :: elem
  type(tem_stencilElement_type) :: tStencil
  type(tem_debug_type) :: debug
  integer :: QQ, elemPos, iQQN, addedPos, nElems, iElem, iSize, iStencilSize
  integer, allocatable :: elemSize(:), stencilSize(:)
  logical :: error, wasAdded
  character(len=64) :: buffer
  integer(kind=long_k) :: nBytes, nReqBytes
  integer :: nElems_fluid
  ! ----------------------------------------------------------------------------
  nElems_fluid = 0
  error = .false.
  ! initialize the debug type
  call tem_logging_init( me    = debug%logger,                                 &
    &                    level = 10,                                           &
    &                    rank  = 0,                                            &
    &                    filename = 'tem_elem_memTrace.res' )

  ! Get initial memory demand
  call tem_reportStatus( debug = debug, level = 1, &
                         text = 'starting overhead', string = 'VmRSS')
  allocate(elemSize(1))
  elemSize = [ 100 ] !, 1000000, 16777216 ]
  allocate(stencilSize(1))
  stencilSize = [ 0 ]

  do iStencilSize = 1, size( stencilSize )
    QQ = stencilSize( iStencilSize )
    write(*,*) 'Starting the memory analysis for QQ = ', QQ

    call init( me = tStencil, QQN = QQ-1, headerPos = 1 )
    do iSize = 1, size(elemSize)
      nElems = elemSize( iSize )
      write(*,*) 'analyzing with nElems ', nElems
      call init( me = elem )
      write(buffer,'(a,i10,a,i3)') 'starting to add nElems ', nElems, ' QQ', QQ
      call tem_reportStatus( debug = debug, level = 1, &
                             text = buffer, string = 'VmRSS')
      call getSize( me = elem, elemSize = nBytes )
      call getReqSize( me = elem, elemSize = nReqBytes )
      !write(*,*) 'elemSize at beginning', nBytes, elem%tID%containerSize, nReqBytes

      ! Create a list of dummy elements
      do iElem = 1, nElems
        !write(*,*) 'element number ', iElem
        ! Append a dummy element
        call append( me = elem, &
          &          tID = int( iElem, kind=long_k),  &
          &          pntTID = iElem,                  &
          &          eType = eT_fluid,                &
          &          nNeighIDs = QQ-1,                &
          &          sourceProc = 1,                  &
          &          property = 0_long_k,             &
          &          pos = elemPos, wasAdded = wasAdded )
        ! Append stencil neighbors to the neighbor array of the element
        do iQQN = 1, QQ-1
          call append( me  = elem%neighID%val( elemPos ),         &
            &          val = int( iQQN, kind=long_k),             &
            &          pos = addedPos, wasAdded = wasAdded )
          tStencil%tIDpos( iQQN ) = addedPos
        end do
        if ( QQ > 0 ) then
          ! Append the temporary stencil to the element
          call append( me = elem%stencil%val( elemPos ),            &
            &          val = tStencil )
        end if
      end do

      write(buffer,'(a,i10,a,i3)') 'finished adding nElems ', nElems, ' QQ', QQ
      call tem_reportStatus( debug = debug, level = 1, &
                             text = buffer, string = 'VmRSS')
      call getSize( me = elem, elemSize = nBytes )
      call getReqSize( me = elem, elemSize = nReqBytes )
      write(*,*) 'elemSize after adding all elems', nBytes/nElems,  &
      nReqBytes/nElems  

      call destroy( me = elem )
      write(buffer,'(a,i10,a,i3)') 'after destruction'
      call tem_reportStatus( debug = debug, level = 1, &
                             text = buffer, string = 'VmRSS')
    end do
    call destroy( me = tStencil )
  end do

  close( debug%logger%funit(0) )
  if( error ) then
    write(*,*) 'FAILED'
    stop -1
  else
    write(*,*) 'PASSED'
    stop 0 
  end if

end program elem_memAlloc
