! Copyright (c) 2011-2013, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2011-2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2011 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2011 Jens Zudrop <j.zudrop@grs-sim.de>
! Copyright (c) 2011-2012 Khaled Ibrahim <k.ibrahim@grs-sim.de>
! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2014 Peter Vitt <peter.vitt2@uni-siegen.de>
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
! ****************************************************************************** !
!> Module to provide smart growing data structures.
!!
!! The dynamic arrays provided by this module are
!! capable of handling lists of values, which might
!! need to grow over time.
!! Removal of entries is not possible directly.
!! A ranking list for the indices is maintained to
!! fast lookups of given values by a binary search.
!! This allows for the efficient handling of lists
!! with unique entries.
!! The complete module might be put into a CoCo Text
!! template, to create new modules of this object
!! for different types. For now, two different
!! templates are used for the declaration part and
!! the implementation part.
!!
?? include 'arrayMacros.inc'
!!
module tem_dyn_array_module

  ! include treelm modules
  use env_module, only: long_k, rk, minLength, labelLen, zeroLength

  implicit none

! -----------------------------------------------------------------
! In order to define this data structure easily for several data
! type, we use the Coco text copying feature here, to duplicate
! the necessary declarations.
! tname: indicates type of dynamic array (long, int, real, ...)
! tstring: is the actual string describing the type specification

! Here the actual declarations are put in by CoCo:
?? copy :: DA_decltxt(long, integer(kind=long_k))
?? copy :: DA_decltxt(int, integer)
?? copy :: DA_decltxt(real, real(kind=rk))
?? copy :: DA_decltxt(label, character(len=labelLen))
! -----------------------------------------------------------------

contains

! -----------------------------------------------------------------
! Also for the implementation, we use the copy feature, to provide
! the necessary duplications to deal with the various types.
! tname ... indicates type of dynamic array (long, int, real, ...)

! ****************************************************************************** !
?? copy :: DA_impltxt(long, integer(kind=long_k), integer(kind=long_k))
?? copy :: DA_impltxt(int, integer, integer)
?? copy :: DA_impltxt(real, real(kind=rk), real(kind=rk))
?? copy :: DA_impltxt(label, character(len=labelLen), character(len=*))
! ****************************************************************************** !


end module tem_dyn_array_module
! ****************************************************************************** !
!> Dynamic data structures
!!
!! During the processing of the mesh to create the solver data structures,
!! common tasks arise.
!! These include among others
!! a unique collection of elements in lists, i.e. without duplicates,
!! a localization of a stored element,
!! sorting of a list of records.
!! As Fortran does not intrinsically include dynamically growing data
!! structures, they first have to be introduced manually.
!!
!! \section tem_dyn_retrievePos Retrieving the position of a record
!! You might want to retrieve the position of a (unique) value within a list.
!! In order to achieve this efficiently, the list with \f$n\f$ entries must be
!! ordered so we can rely on the binary search, which performs with
!! \f$\mathcal{O}( \log n )\f$.
!!
!! \section tem_dyn_dataStructures Dynamic data structures
!! Hence we introduce a derived data type which allows an efficient and simple
!! element adding and retrieving process, as well as a sorted order of the
!! entries.
!! A [module] (@ref tem_dyn_array_module) provides the type definitions and all
!! the related procedures to act on this data type.
!! The data type includes the size of the container within, the amount of its
!! actual entries, a flag for uniqueness and the container array itself as
!! well as a conforming array for the sorted retrieval of the entries in the
!! value array.
!!
!! \todo What shall be listed in lst:tem_dyn_array ??
!!
!! The definition is depicted in listing \ref lst:tem_dyn_array where macros
!! are used for the label of a dynamic type `?tname?`and the used data type
!! `?tstring?`.
!! The introduction of macros allow for a generic
!! definition for all data types of the container,
!! whereas the compiled source code contains the definitions for each required
!! data type.
!!
!! In order to use such data structures they first have to be initialized with
!! options such as the uniqueness of the entries or the initial size of the
!! container.
!! The initialization simply allocates the container with the given options and
!! resets the number of records.
!! Correspondingly, the object can be destroyed by `destroy( me = valName )`.
!!
!!
!! If an object has been initialized, records can be added by simply calling
!! `append`.
!! Subsequently, checks are performed, for the existence of the entry in the
!! list and the sufficiency of the size of the container to incorporate another
!! entry.
!! In case the size does not suffice, the array is doubled in size through a
!! temporary array, and then the record is added.
!! Information about if and where the entry was stored is returned to the
!! caller.
