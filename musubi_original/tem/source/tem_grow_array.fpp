! Copyright (c) 2011-2014 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2014-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
! ****************************************************************************** !
!> Module to provide simple growing data structures.
!! The dynamic arrays provided by this module are
!! capable of handling lists of values, which might
!! need to grow over time.
!! Removal of entries is not possible directly.
!! The complete module might be put into a CoCo Text
!! template, to create new modules of this object
!! for different types. For now, two different
!! templates are used for the declaration part and
!! the implementation part.
?? include 'arrayMacros.inc'
!!
module tem_grow_array_module

  ! include treelm modules
  use env_module, only: long_k, rk, minLength, zeroLength, labelLen

  implicit none

  type intArray2d_type
    integer, allocatable :: val(:,:)
  end type

! -----------------------------------------------------------------
! Growing array (GA)
! tname ... indicates type of dynamic array (long, int, real, ...)

?? copy :: GA_decltxt(logical, logical)
?? copy :: GA_decltxt(long, integer(kind=long_k))
?? copy :: GA_decltxt(int, integer)
?? copy :: GA_decltxt(real, real(kind=rk))
?? copy :: GA_decltxt(dtInt2d, type(intArray2d_type))
?? copy :: GA_decltxt(label, character(len=labelLen))
?? copy :: GA_decltxt(char, character)

! Growing array
! -----------------------------------------------------------------


! -----------------------------------------------------------------
! Growing array (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)
! 2d array, but only the second dimension can grow.
!

?? copy :: GA2d_decltxt(char, character)
?? copy :: GA2d_decltxt(logical, logical)
?? copy :: GA2d_decltxt(long, integer(kind=long_k))
?? copy :: GA2d_decltxt(int, integer)
?? copy :: GA2d_decltxt(real, real(kind=rk))

! Growing array 2D
! -----------------------------------------------------------------




contains

! -----------------------------------------------------------------
! Growing array (GA)
! tname ... indicates type of dynamic array (long, int, real, ...)

?? copy :: GA_impltxt(logical, logical)
?? copy :: GA_impltxt(long, integer(kind=long_k))
?? copy :: GA_impltxt(int, integer)
?? copy :: GA_impltxt(real, real(kind=rk))
?? copy :: GA_impltxt(dtInt2d, type(intArray2d_type))
?? copy :: GA_impltxt(label, character(len=labelLen), character(len=*))
?? copy :: GA_impltxt(char, character)
! Growing array
! -----------------------------------------------------------------



! ****************************************************************************** !
! -----------------------------------------------------------------
! 2d Array, which can grow in second dimension only (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)

?? copy :: GA2d_impltxt(char, character, '')
?? copy :: GA2d_impltxt(logical, logical, .false.)
?? copy :: GA2d_impltxt(long, integer(kind=long_k), -1_long_k)
?? copy :: GA2d_impltxt(int, integer, -1)
?? copy :: GA2d_impltxt(real, real(kind=rk), -1.0_rk)

! Growing array 2D
! -----------------------------------------------------------------


end module tem_grow_array_module
! ****************************************************************************** !
