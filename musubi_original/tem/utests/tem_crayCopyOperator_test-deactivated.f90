! Copyright (c) 2014 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
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
program test_copyOperator
  implicit none
  
  type variable
    integer, allocatable :: val(:) 
  end type variable

  type dyn_variable
    type(variable), allocatable :: var(:)
  end type dyn_variable

  type(variable) :: var1
  type(dyn_variable) :: dynVar

  allocate(dynVar%var(1))
  allocate(dynVar%var(1)%val(1))
  dynVar%var(1)%val = (/ 1 /)
  write(*,*) 'PASSED 1', ' dynVar%var1 size ', size(dynVar%var(1)%val), &
    &        ' allocated ', allocated(dynVar%var(1)%val)
  allocate(var1%val(1))
  var1 = dynVar%var(1)
  write(*,*) 'PASSED 2', ' var1 size ', size(var1%val), &
    &        ' allocated ', allocated(var1%val)

end program test_copyOperator
