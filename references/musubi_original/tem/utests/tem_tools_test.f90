! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
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
program tem_tools_test

  use env_module,       only: rk, long_k
  use tem_tools_module, only: tem_getOptValOrDef
  use tem_float_module, only: operator(.feq.)

  implicit none

  logical :: res = .true.

  res = res .and. test_character() == 'X'
  res = res .and. test_character('Y') =='Y'
  res = res .and. test_int() == 42
  res = res .and. test_int(23) == 23
  res = res .and. test_logical() .eqv. .true.
  res = res .and. test_logical(.false.) .eqv. .false.
  res = res .and. test_long() == 42_long_k
  res = res .and. test_long(23_long_k) == 23_long_k
  res = res .and. (test_real() .feq. 4.2_rk)
  res = res .and. (test_real(2.3_rk) .feq. 2.3_rk)

  if(res) write(*,*) 'PASSED'

contains

  function test_character(char) result(res)
    character :: res
    character, intent(in), optional :: char
    res = tem_getOptValOrDef(char,'X')
  end function

  function test_int(int) result(res)
    integer :: res
    integer, intent(in), optional :: int
    res = tem_getOptValOrDef(int,42)
  end function

  function test_logical(bool) result(res)
    logical :: res
    logical, intent(in), optional :: bool
    res = tem_getOptValOrDef(bool,.true.)
  end function

  function test_long(long) result(res)
    integer(kind=long_k) :: res
    integer(kind=long_k), intent(in), optional :: long
    res = tem_getOptValOrDef(long,42_long_k)
  end function

  function test_real(float) result(res)
    real(kind=rk) :: res
    real(kind=rk), intent(in), optional :: float
    res = tem_getOptValOrDef(float,4.2_rk)
  end function

end program tem_tools_test
