! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2013, 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! This program tries to evaluate the memory demand for derived types
!
program derivedType_memAlloc

  ! include treelm modules
  use tem_debug_module,         only: tem_reportStatus, tem_debug_type
  use tem_dyn_array_module,     only: dyn_longArray_type, init, append, destroy
  use tem_arrayofarrays_module, only: grw_dynlongArray_type, init, append, &
    &                                 destroy
  use tem_logging_module,       only: tem_logging_init

  implicit none
  type innerDT_type
    integer :: test
  end type innerDT_type

  type outerDT_type
    type( innerDT_type ) :: inner
  end type outerDT_type
  type outerDTarr_type
    type( innerDT_type ), allocatable :: inner(:)
  end type outerDTarr_type
  type outerDTpnt_type
    type( innerDT_type ), pointer :: inner(:)
  end type outerDTpnt_type
  !----------------------------------------------------------------------------
  type(outerDT_type), allocatable :: simpleArr(:)
  type(outerDTarr_type), allocatable :: nestedArr(:)
  type(outerDTpnt_type), allocatable :: nestedPnt(:)
  type(grw_dynlongArray_type) :: array
  type(dyn_longArray_type) :: tDynLongArray
  type(tem_debug_type) :: debug
  integer :: QQ, nElems, iElem, iSize
  integer :: iArraySize, performTest, iInner
  integer, allocatable :: elemSize(:), arraySize(:)
  logical :: error
  character(len=64) :: buffer
  !----------------------------------------------------------------------------
  performTest = 4
  error = .false.

  ! initialize the debug type
  call tem_logging_init( me    = debug%logger,                                 &
    &                    level = 10,                                           &
    &                    rank  = 0,                                            &
    &                    filename = 'tem_dT_memTrace.res' )

  ! Get initial memory demand
  call tem_reportStatus( debug = debug, level = 1, &
                         text = 'starting overhead', string = 'VmRSS')
  allocate(elemSize(3))
  elemSize = [ 10000, 100000, 1000000 ]
  allocate(arraySize(1))
  arraySize = [ 10 ]

  select case( performTest )
  case default
    write(*,*) 'error, no valid test chosen.'
    write(*,*) 'choose: 1 ... grw_dynlongArray_test'
    write(*,*) '        2 ... simple derived type test'
    write(*,*) '        3 ... nested derived type test'

  case( 4 ) ! nested derived type test
    do iArraySize = 1, size( arraySize )
      QQ = arraySize( iArraySize )
      write(*,*) 'Starting the memory analysis for internal arr size of ', QQ
      do iSize = 1, size(elemSize)
        nElems = elemSize( iSize )
        write(*,*) 'Starting the memory analysis for nested pointer derived type arrays'
        write(buffer,'(a,i10,a,i3)') 'starting to allocate with nElems ', nElems
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')
        allocate( nestedPnt( nElems ))
        do iElem = 1, nElems
          allocate( nestedPnt( iElem )%inner( QQ ))
          do iInner = 1, QQ
            nestedPnt( iElem )%inner( iInner )%test = 1
          enddo
        enddo
        write(buffer,'(a,i10,a,i3)') 'finished allocation    ', nElems
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')

        do iElem = 1, nElems
          deallocate( nestedPnt( iElem )%inner )
        end do
        deallocate( nestedPnt )
        write(buffer,'(a,i10,a,i3)') 'after deallocation'
        call tem_reportStatus( debug = debug, level = 1, &
                                 text = buffer, string = 'VmRSS')
      end do
    end do
  case( 3 ) ! nested derived type test
    do iArraySize = 1, size( arraySize )
      QQ = arraySize( iArraySize )
      write(*,*) 'Starting the memory analysis for internal arr size of ', QQ
      do iSize = 1, size(elemSize)
        nElems = elemSize( iSize )
        write(*,*) 'Starting the memory analysis for nested allocatable derived type arrays'
        write(buffer,'(a,i10,a,i3)') 'starting to allocate with nElems ', nElems
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')
        allocate( nestedArr( nElems ))
        do iElem = 1, nElems
          allocate( nestedArr( iElem )%inner( QQ ))
          do iInner = 1, QQ
            nestedArr( iElem )%inner( iInner )%test = 1
          enddo
        enddo
        write(buffer,'(a,i10,a,i3)') 'finished allocation    ', nElems
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')

        do iElem = 1, nElems
          deallocate( nestedArr( iElem )%inner )
        end do
        deallocate( nestedArr )
        write(buffer,'(a,i10,a,i3)') 'after deallocation'
        call tem_reportStatus( debug = debug, level = 1, &
                                 text = buffer, string = 'VmRSS')
      end do
    end do
  case( 2 ) ! simple derived type test
    do iSize = 1, size(elemSize)
      nElems = elemSize( iSize )
      write(*,*) 'Starting the memory analysis for simple nested derived type arrays'
      write(buffer,'(a,i10,a,i3)') 'starting to allocate with nElems ', nElems
      call tem_reportStatus( debug = debug, level = 1, &
                             text = buffer, string = 'VmRSS')
      allocate( simpleArr( nElems ))
      do iElem = 1, nElems
        simpleArr( iElem )%inner%test = 1
      enddo
      write(buffer,'(a,i10,a,i3)') 'finished allocation    ', nElems
      call tem_reportStatus( debug = debug, level = 1, &
                             text = buffer, string = 'VmRSS')

      deallocate( simpleArr )
      write(buffer,'(a,i10,a,i3)') 'after deallocation'
      call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')

    end do
  case( 1 ) ! grw_dynlongArray test
    do iArraySize = 1, size( arraySize )
      QQ = arraySize( iArraySize )
      write(*,*) 'Starting the memory analysis for internal arr size of ', QQ

      do iSize = 1, size(elemSize)
        nElems = elemSize( iSize )
        write(*,*) 'analyzing with nElems ', nElems
        call init( me = array )
        write(buffer,'(a,i10,a,i3)') 'starting to add nElems ', nElems, ' QQ', QQ
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')

        ! Create a list of dummy entries
        do iElem = 1, nElems
          !write(*,*) 'element number ', iElem
          ! Append a dummy element
          call append( me = array, val = tDynLongArray )
          array%val( iElem )%nVals = 0
        end do

        write(buffer,'(a,i10,a,i3)') 'finished adding nElems ', nElems, ' QQ', QQ
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')

        call destroy( me = array )
        write(buffer,'(a,i10,a,i3)') 'after destruction'
        call tem_reportStatus( debug = debug, level = 1, &
                               text = buffer, string = 'VmRSS')
      end do
      call destroy( me = tDynLongArray )
    end do
  end select

  close( debug%logger%funit(0) )

  if( error ) then
    write(*,*) 'FAILED'
    stop -1
  else
    write(*,*) 'PASSED'
    stop 0
  end if

end program derivedType_memAlloc
