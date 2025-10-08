! Copyright (c) 2016, 2019 Peter Vitt <peter.vitt2@uni-siegen.de>
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
program tem_float_test

  use env_module,       only: rk
  use tem_float_module, only: operator(.feq.), &
    &                         operator(.fne.), &
    &                         operator(.fgt.), &
    &                         operator(.fge.), &
    &                         operator(.flt.), &
    &                         operator(.fle.)

  implicit none

  real(kind=rk) :: testvalue
  real(kind=rk) :: testarray(2)
  logical :: res = .true.

  write(*,*) 'testvalue = 0._rk'
  testvalue = 0._rk
  testarray = (/ 0._rk, 1._rk /)

  ! 0 == 0
  res = res .and. (testvalue .feq. testvalue)
  call checkTest(res, '(testvalue .feq. testvalue)')
  ! 0 != 0 + something very small
  res = res .and..not. (testvalue .feq. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .feq. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue-spacing(testvalue)))')
  ! 0 + something very small != 0
  res = res .and..not. ((testvalue+spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .feq. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .feq. testvalue)')

  ! (0,1) == (0,1)
  res = res .and. (testarray .feq. testarray)
  call checkTest(res, '(testarray .feq. testarray)')

  ! !0 != 0
  res = res .and..not. (testvalue .fne. testvalue)
  call checkTest(res, '.not. (testvalue .fnq. testvalue)')
  ! !0 == 0 + something very small
  res = res .and. (testvalue .fne. (testvalue+spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue+spacing(testvalue)))')
  res = res .and. (testvalue .fne. (testvalue-spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue-spacing(testvalue)))')
  ! !0 + something very small == 0
  res = res .and. ((testvalue+spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue+spacing(testvalue)) .fne. testvalue)')
  res = res .and. ((testvalue-spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue-spacing(testvalue)) .fne. testvalue)')

  ! (0,1) != (0,1)
  res = res .and..not. (testarray .fne. testarray)
  call checkTest(res, '.not. (testarray .fne. testarray)')

  ! 0 lt 0 + something very small
  res = res .and. (testvalue .flt. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .flt. (testvalue+spacing(testvalue)))')
  ! 0 !gt 0 + something very small
  res = res .and..not. (testvalue .fgt. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fgt. (testvalue+spacing(testvalue)))')
  ! 0 + something very small !lt 0
  res = res .and..not. ((testvalue+spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .flt. testvalue)')
  ! 0 + something very small gt 0
  res = res .and. ((testvalue+spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fgt. testvalue)')
  ! 0 !lt 0 - something very small
  res = res .and..not. (testvalue .flt. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .flt. (testvalue-spacing(testvalue)))')
  ! 0 gt 0 - something very small
  res = res .and. (testvalue .fgt. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fgt. (testvalue-spacing(testvalue)))')
  ! 0 - something very small lt 0
  res = res .and. ((testvalue-spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .flt. testvalue)')
  ! 0 - something very small !gt 0
  res = res .and..not. ((testvalue-spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fgt. testvalue)')

  ! 0 lte 0
  res = res .and. (testvalue .fle. testvalue)
  call checkTest(res, '(testvalue .fle. testvalue)')
  ! 0 lge 0
  res = res .and. (testvalue .fge. testvalue)
  call checkTest(res, '(testvalue .fge. testvalue)')
  ! 0 lte 0 + something very small
  res = res .and. (testvalue .fle. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .fle. (testvalue+spacing(testvalue)))')
  ! 0 !gte 0 + something very small
  res = res .and..not. (testvalue .fge. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fge. (testvalue+spacing(testvalue)))')
  ! 0 + something very small !lte 0
  res = res .and..not. ((testvalue+spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .fle. testvalue)')
  ! 0 + something very small gte 0
  res = res .and. ((testvalue+spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fge. testvalue)')
  ! 0 !lte 0 - something very small
  res = res .and..not. (testvalue .fle. (testvalue-spacing(testvalue)))
  call checkTest(res, '..not. (testvalue .fle. (testvalue-spacing(testvalue)))')
  ! 0 gte 0 - something very small
  res = res .and. (testvalue .fge. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fge. (testvalue-spacing(testvalue)))')
  ! 0 - something very small lte 0
  res = res .and. ((testvalue-spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .fle. testvalue)')
  ! 0 - something very small !gte 0
  res = res .and..not. ((testvalue-spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fge. testvalue)')

  ! The same with huge
  write(*,*) 'testvalue = huge(0._rk)'
  testvalue = huge(0._rk)
  res = res .and. (testvalue .feq. testvalue)
  call checkTest(res, '(testvalue .feq. testvalue)')
  res = res .and..not. (testvalue .feq. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue-spacing(testvalue)))')
  res = res .and..not. ((testvalue-spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .feq. testvalue)')

  res = res .and..not. (testvalue .fne. testvalue)
  call checkTest(res, '.not. (testvalue .fnq. testvalue)')
  res = res .and. (testvalue .fne. (testvalue-spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue-spacing(testvalue)) .fne. testvalue)')

  res = res .and..not. (testvalue .flt. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .flt. (testvalue-spacing(testvalue)))')
  res = res .and. (testvalue .fgt. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fgt. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .flt. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fgt. testvalue)')

  res = res .and. (testvalue .fle. testvalue)
  call checkTest(res, '(testvalue .fle. testvalue)')
  res = res .and. (testvalue .fge. testvalue)
  call checkTest(res, '(testvalue .fge. testvalue)')
  res = res .and..not. (testvalue .fle. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fle. (testvalue-spacing(testvalue)))')
  res = res .and. (testvalue .fge. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fge. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .fle. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fge. testvalue)')

  ! and now with tiny
  write(*,*) 'testvalue = tiny(0._rk)'
  testvalue = tiny(0._rk)
  res = res .and. (testvalue .feq. testvalue)
  call checkTest(res, '(testvalue .feq. testvalue)')
  res = res .and..not. (testvalue .feq. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .feq. testvalue)')

  res = res .and..not. (testvalue .fne. testvalue)
  call checkTest(res, '.not. (testvalue .fnq. testvalue)')
  res = res .and. (testvalue .fne. (testvalue+spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue+spacing(testvalue)))')
  res = res .and. (testvalue .fne. (testvalue-spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue+spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue+spacing(testvalue)) .fne. testvalue)')
  res = res .and. ((testvalue-spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue-spacing(testvalue)) .fne. testvalue)')

  res = res .and. (testvalue .flt. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .flt. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .fgt. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fgt. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .flt. testvalue)')
  res = res .and. ((testvalue+spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fgt. testvalue)')

  res = res .and. (testvalue .fle. testvalue)
  call checkTest(res, '(testvalue .fle. testvalue)')
  res = res .and. (testvalue .fge. testvalue)
  call checkTest(res, '(testvalue .fge. testvalue)')
  res = res .and. (testvalue .fle. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .fle. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .fge. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fge. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .fle. testvalue)')
  res = res .and. ((testvalue+spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fge. testvalue)')

  ! Now with some rather big value
  write(*,*) 'testvalue = 1234567890_rk * 10**9'
  testvalue = 1234567890_rk * 10**9
  res = res .and. (testvalue .feq. testvalue)
  call checkTest(res, '(testvalue .feq. testvalue)')
  res = res .and..not. (testvalue .feq. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .feq. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue-spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .feq. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .feq. testvalue)')

  res = res .and..not. (testvalue .fne. testvalue)
  call checkTest(res, '.not. (testvalue .fnq. testvalue)')
  res = res .and. (testvalue .fne. (testvalue+spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue+spacing(testvalue)))')
  res = res .and. (testvalue .fne. (testvalue-spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue+spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue+spacing(testvalue)) .fne. testvalue)')
  res = res .and. ((testvalue-spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue-spacing(testvalue)) .fne. testvalue)')

  res = res .and. (testvalue .flt. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .flt. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .fgt. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fgt. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .flt. testvalue)')
  res = res .and. ((testvalue+spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fgt. testvalue)')
  res = res .and..not. (testvalue .flt. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .flt. (testvalue-spacing(testvalue)))')
  res = res .and. (testvalue .fgt. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fgt. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .flt. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fgt. testvalue)')

  res = res .and. (testvalue .fle. testvalue)
  call checkTest(res, '(testvalue .fle. testvalue)')
  res = res .and. (testvalue .fge. testvalue)
  call checkTest(res, '(testvalue .fge. testvalue)')
  res = res .and. (testvalue .fle. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .fle. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .fge. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fge. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .fle. testvalue)')
  res = res .and. ((testvalue+spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fge. testvalue)')
  res = res .and..not. (testvalue .fle. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fle. (testvalue-spacing(testvalue)))')
  res = res .and. (testvalue .fge. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fge. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .fle. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fge. testvalue)')

  ! and with some rather small value
  write(*,*) 'testvalue = 9876543210_rk * 10_rk**(-100)'
  testvalue = 9876543210._rk * 10._rk**(-100)
  res = res .and. (testvalue .feq. testvalue)
  call checkTest(res, '(testvalue .feq. testvalue)')
  res = res .and..not. (testvalue .feq. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .feq. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .feq. (testvalue-spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .feq. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .feq. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .feq. testvalue)')

  res = res .and..not. (testvalue .fne. testvalue)
  call checkTest(res, '.not. (testvalue .fnq. testvalue)')
  res = res .and. (testvalue .fne. (testvalue+spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue+spacing(testvalue)))')
  res = res .and. (testvalue .fne. (testvalue-spacing(testvalue)))
  call checkTest(res, ' (testvalue .fne. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue+spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue+spacing(testvalue)) .fne. testvalue)')
  res = res .and. ((testvalue-spacing(testvalue)) .fne. testvalue)
  call checkTest(res, ' ((testvalue-spacing(testvalue)) .fne. testvalue)')

  res = res .and. (testvalue .flt. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .flt. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .fgt. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fgt. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .flt. testvalue)')
  res = res .and. ((testvalue+spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fgt. testvalue)')
  res = res .and..not. (testvalue .flt. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .flt. (testvalue-spacing(testvalue)))')
  res = res .and. (testvalue .fgt. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fgt. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .flt. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .flt. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .fgt. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fgt. testvalue)')

  res = res .and. (testvalue .fle. testvalue)
  call checkTest(res, '(testvalue .fle. testvalue)')
  res = res .and. (testvalue .fge. testvalue)
  call checkTest(res, '(testvalue .fge. testvalue)')
  res = res .and. (testvalue .fle. (testvalue+spacing(testvalue)))
  call checkTest(res, '(testvalue .fle. (testvalue+spacing(testvalue)))')
  res = res .and..not. (testvalue .fge. (testvalue+spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fge. (testvalue+spacing(testvalue)))')
  res = res .and..not. ((testvalue+spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '.not. ((testvalue+spacing(testvalue)) .fle. testvalue)')
  res = res .and. ((testvalue+spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '((testvalue+spacing(testvalue)) .fge. testvalue)')
  res = res .and..not. (testvalue .fle. (testvalue-spacing(testvalue)))
  call checkTest(res, '.not. (testvalue .fle. (testvalue-spacing(testvalue)))')
  res = res .and. (testvalue .fge. (testvalue-spacing(testvalue)))
  call checkTest(res, '(testvalue .fge. (testvalue-spacing(testvalue)))')
  res = res .and. ((testvalue-spacing(testvalue)) .fle. testvalue)
  call checkTest(res, '((testvalue-spacing(testvalue)) .fle. testvalue)')
  res = res .and..not. ((testvalue-spacing(testvalue)) .fge. testvalue)
  call checkTest(res, '.not. ((testvalue-spacing(testvalue)) .fge. testvalue)')

  if ( res ) write(*,*) 'PASSED'

contains

  subroutine checkTest(res, msg)
    logical, intent(in) :: res
    character(len=*), intent(in) :: msg

    if ( .not. res) then
      write(*,*) 'Test failed at ' // msg
      stop
    end if
  end subroutine

end program tem_float_test
