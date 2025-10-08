! Copyright (c) 2012-2013, 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2014, 2017 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2014 Simon Zimny <s.zimny@grs-sim.de>
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
program bin_search_test
  use env_module, only: long_k, rk
  use tem_timer_module, only: tem_addTimer, tem_startTimer, tem_stopTimer,     &
    &                         tem_labeledtimer_type, tem_getTimerVal
  use tem_dyn_array_module, only: dyn_longArray_type, init, append,            &
    &                             PositionOfVal
  use tem_general_module,    only: tem_general_type, tem_start, tem_finalize

  !mpi!nprocs = 1

  implicit none

  integer, parameter :: nVals = 1000000

  integer :: long_timing
  integer :: iVal
  integer, allocatable :: pos(:)
  type(dyn_longArray_type) :: longArray
  type(tem_labeledtimer_type) :: timer
  type(tem_general_type) :: general
  real(kind=rk) :: timerVal

  allocate(pos(nVals))

  write(*,*) 'Running binary search test...'
  write(*,*) 'on ', nVals, ' Values'

  ! Init the Treelm environment
  call tem_start('TREELM unit test', general)

  call tem_addTimer( me          = timer,        &
    &                timerHandle = long_timing,  &
    &                timerName   = 'long_timing' )

  call init(me = longArray, length = nVals)

  call tem_startTimer( timer%timedat, long_timing )
  do iVal=1,nVals
    call append(me = longArray, val = int(iVal, kind=long_k), pos = pos(iVal))
  end do
  call tem_stopTimer( timer%timedat, long_timing )
  timerVal = tem_getTimerVal(timer%timedat, long_timing)
  write(*,'(a,f15.4,a)') ' Run-time long_timing:', timerVal, ' s'

  call tem_finalize(general)

  ! access an element from the pos array, to stop compiler from optimizing the
  ! loop away:
  write(*,*) 'last Value has position:', pos(nVals)

  write(*,*) 'PASSED'

end program bin_search_test
