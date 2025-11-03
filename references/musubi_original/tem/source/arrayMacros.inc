! Copyright (c) 2012-2013 Simon Zimny <s.zimny@grs-sim.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2012-2016 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012, 2015-2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2013 Daniel Harlacher <d.harlacher@grs-sim.de>
! Copyright (c) 2014 Kartik Jain <kartik.jain@uni-siegen.de>
! Copyright (c) 2014 Verena Krupp <verena.krupp@uni-siegen.de>
! Copyright (c) 2015-2017 Jiaxing Qi <jiaxing.qi@uni-siegen.de>
! Copyright (c) 2015-2016 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Daniel Fleischer <daniel.fleischer@student.uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2017 Daniel Petró <daniel.petro@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach, Simon Zimny and Manuel
! Hasert for German Research School for Simulation Sciences GmbH.
!
! Parts of this file were written by Harald Klimach, Kannan Masilamani,
! Daniel Harlacher, Kartik Jain, Verena Krupp, Jiaxing Qi, Peter Vitt,
! Daniel Fleischer, Tobias Girresser and Daniel Petró for University Siegen.
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

! This file contains the source code for growing and dynamic arrays.
! This is used for arrays of primitives (int, long_int, real, ...) as well as
! for arrays of derived datatypes (tem_variable_type,...).
!
! To use these macros include the following to your source file.
!
! Smart growing array (GA) for ?tstring?
! Growing Arrays:
!
! declaration
!
?? text :: GA_decltxt(tname, tstring)
  !> Growing array type for ?tstring?
  type grw_?tname?Array_type
    integer :: nVals = 0
    integer :: ContainerSize = 0
    ?tstring?, allocatable :: Val(:)
  end type

  !> Initialize the dynamic array
  interface init
    module procedure init_GA_?tname?
  end interface

  !> Truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_GA_?tname?
  end interface

  !> Empty the entries  without changing arrays
  interface empty
    module procedure empty_GA_?tname?
  end interface

  !> Destroy the dynamic array
  interface destroy
    module procedure destroy_GA_?tname?
  end interface

  !> Insert an element at a given position
  interface placeAt
    module procedure placeAt_GA_?tname?
    module procedure placeAt_GA_?tname?_vec
  end interface

  !> Append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_GA_?tname?
    module procedure append_GA_?tname?_vec
  end interface

  !> Increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_GA_?tname?
  end interface

?? end text GA_decltxt
!
! implementation
!
?? text :: GA_impltxt(tname, tstring)
  subroutine init_GA_?tname?(me, length)
    type(grw_?tname?Array_type), intent(out) :: me !< Dynamic Array to init
    integer, intent(in), optional :: length !< Initial length of the container

    if (present(length)) then
      me%ContainerSize = length
    else
      me%ContainerSize = zeroLength
    end if
    ! Deallocate ...
    if( allocated( me%Val ))     &
      deallocate(me%Val)
    ! ... and reallocate
    allocate(me%Val(me%ContainerSize))
    me%nVals = 0

  end subroutine init_GA_?tname?

  subroutine destroy_GA_?tname?(me)
    type(grw_?tname?Array_type), intent(inout) :: me !< Dynamic Array to destroy

    me%ContainerSize = 0
    me%nVals = 0
    if( allocated( me%Val ) ) deallocate(me%Val)

  end subroutine destroy_GA_?tname?


  subroutine truncate_GA_?tname?(me)
    !------------------------------------------------------------------------
    type(grw_?tname?Array_type) :: me !< Array to truncate
    !------------------------------------------------------------------------
    ?tstring?, allocatable :: tArray(:)
    !------------------------------------------------------------------------
    integer :: ii
    !------------------------------------------------------------------------

    ! Nothing to do if container size is not larger than the number of values
    ! in the array.
    if (me%containerSize > me%nVals) then
      allocate(tArray(me%nVals))
      do ii = 1, me%nVals
        tArray(ii) = me%val(ii)
      end do
      call move_alloc(tarray, me%val)
      me%containerSize = me%nVals
    end if

  end subroutine truncate_GA_?tname?


  subroutine empty_GA_?tname?(me)
    !------------------------------------------------------------------------
    type(grw_?tname?Array_type) :: me !< Array to sortTruncate
    !------------------------------------------------------------------------

    me%nVals = 0

  end subroutine empty_GA_?tname?

  !> Adds the value to a given position inside the growing array.
  !!
  !! If the requested position is outside the current array bounds, the array
  !! will be resized accordingly. If it is inside the current array bounds, the
  !! element at the requested position will be replaced.
  subroutine placeAt_GA_?tname?(me, val, pos, length)
    type(grw_?tname?Array_type) :: me !< Array to place the value into
    ?tstring?, intent(in) :: val !< Value to place at the given position
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length


    ! Value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! In this case FoundPos = me%nVals + 1 holds.
    if (Pos > me%ContainerSize) then
      ! Expand the array, if its boundary is reached
      call expand(me = me, pos = pos, length = length)
    end if

    me%nVals = max( pos, me%nVals )
    me%val(pos) = val

  end subroutine placeAt_GA_?tname?


  !> Adds the values starting from a given position inside the growing array.
  !!
  !! If the requested position is outside the current array bounds, the array
  !! will be resized accordingly. If it is inside the current array bounds, the
  !! elements starting from the requested position will be replaced up to
  !! the element at position `pos + size(val) - 1`.
  subroutine placeAt_GA_?tname?_vec(me, val, pos, length)
    type(grw_?tname?Array_type) :: me !< Array to append the value to
    ?tstring?, intent(in) :: val(:) !< Values to append
    integer, intent(in) :: pos !< predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: ub, ii

    if (me%nVals == huge(me%nVals)) then
      write(*,*) "Reached end of integer range for growing array!"
      write(*,*) "Aborting!!"
      STOP
    end if

    ub = pos + size(val) - 1

    if (ub > me%ContainerSize) then
      ! Expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nVals = max( ub, me%nVals )
    do ii = pos, ub
      me%val(ii) = val(1+ii-pos)
    end do

  end subroutine placeAt_GA_?tname?_vec


  subroutine append_GA_?tname?(me, val, length)
    type(grw_?tname?Array_type) :: me !< Array to append the value to
    ?tstring?, intent(in) :: val !< Value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    ! Value to append is larger than all existing ones,
    ! just put it to the end of the list, this captures
    ! also the case of empty lists.
    ! In this case FoundPos = me%nVals + 1 holds.
    if (me%nVals+1 > me%ContainerSize) then
      ! Expand the array, if its boundary is reached
      call expand(me = me, length = length)
    end if

    me%nVals = me%nVals+1
    me%val(me%nVals) = val

  end subroutine append_GA_?tname?

  subroutine append_GA_?tname?_vec(me, val, length)
    type(grw_?tname?Array_type) :: me !< Array to append the value to
    ?tstring?, intent(in) :: val(:) !< Values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length

    integer :: lb, ub, ii

    if (me%nVals == huge(me%nVals)) then
      write(*,*) "Reached end of integer range for growing array!"
      write(*,*) "Aborting!!"
      STOP
    end if

    lb = me%nVals + 1
    ub = lb + size(val) - 1

    if (ub > me%ContainerSize) then
      ! Expand the array, if its boundary is reached
      call expand(me = me, pos = ub, length = length)
    end if

    me%nVals = max( ub, me%nVals )
    do ii = lb, ub
      me%val(ii) = val(1+ii-lb)
    end do

  end subroutine append_GA_?tname?_vec


  subroutine expand_GA_?tname?(me, pos, length)
    type(grw_?tname?Array_type) :: me !< Array to resize
    integer, intent(in), optional :: pos !< optional predefined position
    !> optional length to expand the array
    integer, intent(in), optional :: length

    ?tstring?, allocatable :: swpval(:)
    integer :: expLen, ii

    expLen = 0
    ! Increase the container by the requested length of double it
    if( present(length) ) then
      expLen = max( length, minLength )
    else
      ! Set the global minimum length, if doubling would be smaller than that
      expLen = max(me%ContainerSize, minLength)
    end if

    ! If a position is given, increase the container to at least the size to
    ! fit the position.
    if( present(pos) ) expLen = max(expLen, pos-me%ContainerSize)

    ! If the current size plus expLen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%ContainerSize) - expLen) <= me%ContainerSize) then
      ! set max container size
      me%ContainerSize = huge(me%ContainerSize)
    else
      ! Set the new container size
      me%ContainerSize = me%ContainerSize + expLen
    end if

    if ( me%nVals > 0 ) then
      allocate(swpval(me%ContainerSize))
      do ii = 1, me%nVals
        swpval(ii) = me%Val(ii)
      end do
      call move_alloc( swpval, me%Val )
    else ! me%nVals == 0
      if ( allocated(me%Val) ) deallocate( me%Val )
      allocate( me%Val(me%containerSize) )
    end if

  end subroutine expand_GA_?tname?

?? end text GA_impltxt

! -----------------------------------------------------------------
! 2d Array, which can grow in second dimension only (GA2d)
! tname ... indicates type of dynamic array (long, int, real, ...)
?? text :: GA2d_decltxt(tname, tstring)
  type grw_?tname?2dArray_type
    integer :: nVals = 0
    integer :: ContainerSize = 0
    integer :: ContainerWidth = 0
    ?tstring?, allocatable :: Val(:,:)
  end type

  !> Initialize the dynamic array
  interface init
    module procedure init_GA2d_?tname?
  end interface

  !> Destroy the dynamic array
  interface destroy
    module procedure destroy_GA2d_?tname?
  end interface

  !> Append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_SingleGA2d_?tname?
    module procedure append_ArrayGA2d_?tname?
  end interface

  interface placeAt
    module procedure placeAt_SingleGA2d_?tname?
    module procedure placeAt_ArrayGA2d_?tname?
  end interface

  !> Increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_GA2d_?tname?
  end interface

?? end text GA2d_decltxt

?? text :: GA2d_impltxt(tname, tstring, invalidval)
  ! ****************************************************************************
  subroutine init_GA2d_?tname?(me, width, length)
    ! --------------------------------------------------------------------------
    !> Dynamic Array to init
    type(grw_?tname?2dArray_type), intent(out) :: me
    !> width of the container
    integer, intent(in)           :: width
    !> Initial length of the container
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------

    me%ContainerWidth = width
    if (present(length)) me%ContainerSize = length

    allocate(me%Val(me%ContainerWidth,me%ContainerSize))
    ! Reset all values (only positive values are valid)
    me%val(:,:) = ?invalidval?
    me%nVals = 0

  end subroutine init_GA2d_?tname?
  ! ****************************************************************************


  ! ****************************************************************************
  !> Destroy the 2d growing array
  !!
  subroutine destroy_GA2d_?tname?(me)
    ! --------------------------------------------------------------------------
    !> Dynamic Array to init
    type(grw_?tname?2dArray_type), intent(inout) :: me
    ! --------------------------------------------------------------------------

    deallocate(me%val)
    me%nvals  = 0
    me%ContainerSize = 0

  end subroutine destroy_GA2d_?tname?
  ! ****************************************************************************


  ! ****************************************************************************
  !> Append a single value to the growing 2D array.
  !!
  subroutine append_SingleGA2d_?tname?(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> Array to append the value to
    type(grw_?tname?2dArray_type) :: me
    !> Value to append
    ?tstring?, intent(in) :: val
    !> Position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> The position in second dimension the element were added to
    integer, intent(out), optional :: pos2
    ! --------------------------------------------------------------------------
    integer :: newPos
    ! --------------------------------------------------------------------------

    newPos = me%nVals + 1
    if (newPos > me%ContainerSize) then
      ! Expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%val(pos1, newPos) = val

    if( present( pos2 ) ) pos2 = newPos

  end subroutine append_SingleGA2d_?tname?
  ! ****************************************************************************


  ! ****************************************************************************
  !> Append an array of values to the growing 2D array.
  !!
  subroutine append_ArrayGA2d_?tname?(me, val, length, pos)
    ! --------------------------------------------------------------------------
    !> Array to append the value to
    type(grw_?tname?2dArray_type) :: me
    !> Array of values to append
    ?tstring?, intent(in) :: val(:)
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> The position in second dimension the elements were added to
    integer, intent(out), optional :: pos
    ! --------------------------------------------------------------------------
    integer :: newPos
    ! --------------------------------------------------------------------------

    if (me%nVals == huge(me%nVals)) then
       write(*,*) "Reached end of integer range for growing array!"
       write(*,*) "Aborting!!"
       STOP
    end if

    newPos = me%nVals + 1
    if (newPos > me%ContainerSize) then
      ! Expand the array, if its boundary is reached
      call expand( me = me, length = length )
    end if

    me%nVals = max(newPos, me%nVals)
    me%val( : , newPos ) = val(:)

    if( present( pos ) ) pos = newPos

  end subroutine append_ArrayGA2d_?tname?
  ! ****************************************************************************


  ! ****************************************************************************
  !> Append a single value to the growing 2D array.
  !!
  subroutine placeAt_SingleGA2d_?tname?(me, val, pos1, pos2, length)
    ! --------------------------------------------------------------------------
    !> Array to append the value to
    type(grw_?tname?2dArray_type) :: me
    !> Value to append
    ?tstring?, intent(in) :: val
    !> Position in first dimension (cannot grow)
    integer, intent(in) :: pos1
    !> Position in Second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandLength
    ! --------------------------------------------------------------------------

    if (pos2 > me%ContainerSize) then
      expandLength = pos2 - me%ContainerSize
      if( present( length ) ) expandLength = max(expandLength, length)
      ! Expand the array, if its boundary is reached
      call expand( me = me, length = expandLength )
    end if
    me%nVals = max(pos2, me%nVals)
    me%val(pos1,pos2) = val

  end subroutine placeAt_SingleGA2d_?tname?
  ! ****************************************************************************


  ! ****************************************************************************
  !> Append an array of values to the growing 2D array.
  !!
  subroutine placeAt_ArrayGA2d_?tname?(me, val, pos2, length)
    ! --------------------------------------------------------------------------
    !> Array to append the value to
    type(grw_?tname?2dArray_type) :: me
    !> Array of values to append
    ?tstring?, intent(in) :: val(:)
    !> Position in Second dimension (can grow)
    integer, intent(in) :: pos2
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    integer :: expandLength
    ! --------------------------------------------------------------------------

    if (me%nVals == huge(me%nVals)) then
       write(*,*) "Reached end of integer range for growing array!"
       write(*,*) "Aborting!!"
       STOP
    end if
    if (pos2 > me%ContainerSize) then
      expandLength = pos2 - me%ContainerSize
      if( present( length ) ) expandLength = max(expandLength, length)
      ! Expand the array, if its boundary is reached
      call expand( me = me, length = expandLength )
    end if

    me%nVals = max(pos2, me%nVals)
    me%val( : , pos2 ) = val(:)

  end subroutine placeAt_ArrayGA2d_?tname?
  ! ****************************************************************************


  ! ****************************************************************************
  !> Expand the growing 2D array
  subroutine expand_GA2d_?tname?( me, length )
    ! --------------------------------------------------------------------------
    !> Array to resize
    type(grw_?tname?2dArray_type) :: me
    !> optional length to expand the array
    integer, intent(in), optional :: length
    ! --------------------------------------------------------------------------
    ?tstring?, allocatable :: swpval(:,:)
    integer :: expLen
    ! --------------------------------------------------------------------------

    expLen = 0
    ! Increase the container by the requested length of double it
    if( present(length) ) then
      expLen = length
    else
      ! Set the global minimum length, if doubling would be smaller than that
      expLen = max(me%ContainerSize, minLength)
    end if


    ! If the current size plus expLen exceeds the max container size,
    ! reduce the size to the max container size.
    if( (huge(me%ContainerSize) - expLen) <= me%ContainerSize) then
      ! set max container size
      me%ContainerSize = huge(me%ContainerSize)
    else
      ! Set the calculated size
      me%ContainerSize = me%ContainerSize + expLen
    end if


    if ( me%nVals > 0 ) then
      allocate(swpval(me%ContainerWidth,me%ContainerSize))
      !> First reset all entries
      swpval(:,:) = ?invalidval?
      swpval(:,:me%nVals) = me%Val(:,:me%nVals)
      call move_alloc( swpval, me%Val )
    else ! me%nVals == 0
      if ( allocated(me%Val) ) deallocate( me%Val )
      allocate( me%Val(me%ContainerWidth,me%ContainerSize) )
    end if

  end subroutine expand_GA2d_?tname?
  ! ****************************************************************************

! Growing array 2D
! -----------------------------------------------------------------


?? end text GA2d_impltxt

!
!------------------------------------------------------------------------------
!
! dynamic Arrays:
!
! declaration
!
?? text :: DA_decltxt(tname, tstring)
! \brief Smart Dynamic Array (DA) for ?tstring?
!
! This datatype implements a dynamic array, which is capable of
! growing and adding of unique elements. It is available for
! various types, here we deal with $tstring$.
! Sorted array contains the pointers to Val array, instead of
! the actual values in Val array. For example:
! Val:     8, 6, 7, 5
! sorted:  4, 2, 3, 1
  !> dynamic array (DA) type for ?tstring?
  type dyn_?tname?Array_type
    integer :: nVals = 0
    integer :: ContainerSize = 0
    ?tstring?, allocatable :: Val(:)
    integer, allocatable :: sorted(:) !< pointers, NOT values
  end type

  !> Initialize the dynamic array
  interface init
    module procedure init_DA_?tname?
  end interface

  !> Destroy the dynamic array
  interface destroy
    module procedure destroy_DA_?tname?
  end interface

  !> Append a value to the dynamic array
  !! and return its position.
  interface append
    module procedure append_DA_?tname?
    module procedure append_DA_vec?tname?
  end interface

  !> Truncate the array, meaning
  !! cut off the trailing empty entries
  interface truncate
    module procedure truncate_DA_?tname?
  end interface

  !> Empty the array, reset nVals to be 0
  interface empty
    module procedure empty_DA_?tname?
  end interface


  !> Fix the dynamic array, meaning:
  !! store the array in the sorted order and cut
  !! off the trailing empty entries
  interface sortTruncate
    module procedure sortTruncate_DA_?tname?
  end interface

  !> Increase the size of the container
  !! for the array.
  interface expand
    module procedure expand_DA_?tname?
  end interface

  !> Return the position of a given value
  !! in the array val, which is what you usually want to know.
  !! It is the index of a given value
  interface PositionOfVal
    module procedure PosOfVal_?tname?
  end interface

  !> Return the position of a given value
  !! in the list 'sorted'. This is mainly for internal usage.
  !! The sorted list is only a pointer list to the actual values
  !! Thus, in order to get the index of a given value, you
  !! need to look up the entry in the sorted list.
  !! This is done by the PositionOfVal routine
  interface SortedPosOfVal
    module procedure SortPosOfVal_?tname?
  end interface
?? end text DA_decltxt
!
! implementation
!
?? text :: DA_impltxt(tname, tstring, argstring)
! ***************************************************************************** !
  !> Initialization of a dynamic array
  !!
  !! Before a dynamic array can be used, it has to be initialized
  !! with this routine. The initial length provided here, can
  !! avoid reallocations and memory copying, if approximated
  !! correctly enough. If none is specified, the provided container
  !! initially will be of size 0.
  subroutine init_DA_?tname?(me, length)
    !-----------------------------------------------------------------
    type(dyn_?tname?Array_type), intent(out) :: me !< Dynamic Array to init
    integer, intent(in), optional :: length !< Initial length of the container
    !-----------------------------------------------------------------

    if (present(length)) then
      me%ContainerSize = length
    else
      me%ContainerSize = zeroLength
    end if

    ! Deallocate ...
    if( allocated( me%Val ) ) deallocate(me%Val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
    ! ... and reallocate
    allocate(me%Val(me%ContainerSize))
    allocate(me%sorted(me%ContainerSize))
    me%nVals = 0

  end subroutine init_DA_?tname?

  !> Destruction of a dynamic array
  !!
  !! This subroutine takes care of a proper destruction of a
  !! dynamic array, it frees the allocated memory and resets
  !! the internal counts to 0.
  subroutine destroy_DA_?tname?(me)
    type(dyn_?tname?Array_type), intent(inout) :: me !< Dynamic Array to init

    me%ContainerSize = 0
    me%nVals         = 0

    if( allocated( me%Val ) ) deallocate(me%Val)
    if( allocated( me%sorted ) ) deallocate(me%sorted)
  end subroutine destroy_DA_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Appending a value to the dynamic array
  !!
  !! With this subroutine, a given value can be added to the
  !! dynamic array. The actual position of this value in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. With the wasAdded flag, it is indicated,\n
  !! wasAdded = true,  if this entry had to be added,\n
  !! wasAdded = false, if this was already found in the array.
  subroutine append_DA_?tname?(me, val, length, pos, wasAdded )
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: me   !< Array to append the value to
    ?argstring?, intent(in)     :: val  !< Value to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> Position in the array, if the value is found
    integer, intent(out), optional :: pos
    !> Flag to indicate, if val was newly added
    logical, intent(out), optional :: wasAdded
    !------------------------------------------------------------------------
    integer :: FoundPos
    integer :: i
    !------------------------------------------------------------------------

    ! Do a binary search on existing entries (returns closest entry next to
    ! it if not found).
    FoundPos = SortedPosOfVal(me, val, .true.)
    if( present( wasAdded ) ) wasAdded = .false.

    ! If it found the value, the position is smaller than nVals
    if (FoundPos <= me%nVals) then

      ! The returned position might actually be the right entry already or
      ! not, check for this here.
      if ( me%val(me%sorted(FoundPos)) == val ) then

        ! Found the value in a list of unique values,
        ! nothing to do, just return its position.
        if( present( pos ) ) pos = me%sorted(FoundPos)

      else

        ! Need to append a new value!

        if (me%nVals == huge(me%nVals)) then
           write(*,*) "Reached end of integer range for dynamic array!"
           write(*,*) "Aborting!!"
           STOP
        end if

        if( present( wasAdded ) ) wasAdded = .true.
        if (me%nVals == me%ContainerSize) then

          ! Container is full, need to expand it
          call expand(me = me, length = length)
        end if
        me%nVals = me%nVals + 1

        ! Put the new value into the last position in the
        ! array.
        me%val(me%nVals) = val
        do while( FoundPos < me%nVals )
          if(me%val(me%sorted(FoundPos)) /= val) then
            exit
          end if
          ! In case of multiple entries with the same value
          ! move on to the first differing entry.
          FoundPos = FoundPos + 1
        end do
        ! Shift the sorted list of indices, to create a
        ! whole for the value to be inserted, at position
        ! FoundPos.
        do i=me%nVals-1,FoundPos,-1
          me%sorted(i+1) = me%sorted(i)
        end do
        ! Put the index of the new value into the
        ! sorted list at the now freed position.
        me%sorted(FoundPos) = me%nVals
        if( present( pos ) ) pos = me%nVals

      end if

    else

      ! Value to append is larger than all existing ones,
      ! just put it to the end of the list, this captures
      ! also the case of empty lists.
      ! In this case FoundPos = me%nVals + 1 holds.
      if( present( wasAdded ) ) wasAdded = .true.
      if (FoundPos > me%ContainerSize) then
        ! Expand the array, if its boundary is reached
        call expand(me = me, length = length)
      end if
      me%nVals = FoundPos
      me%val(FoundPos) = val
      me%sorted(FoundPos) = FoundPos
      if( present( pos ) ) pos = FoundPos

    end if

  end subroutine append_DA_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Appending a sorted list of values to the dynamic array
  !!
  !! With this subroutine, a given list of sorted values can be added to the
  !! dynamic array. The actual positions of these values in the
  !! dynamic array will be returned, so it can be found again
  !! easily later. With the wasAdded flag, it is indicated,\n
  !! wasAdded = true,  if this entry had to be added,\n
  !! wasAdded = false, if this was already found in the array.
  subroutine append_DA_vec?tname?(me, val, length, pos, wasAdded )
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: me   !< Array to append the value to
    ?argstring?, intent(in)     :: val(:)  !< Values to append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !> Position in the array, the values are found at.
    integer, intent(out), optional :: pos(:)
    !> Flag to indicate, if val was newly added
    logical, intent(out), optional :: wasAdded(:)
    !------------------------------------------------------------------------
    ?tstring? :: lastval
    logical :: addedval(size(val))
    integer :: i
    integer :: veclen
    integer :: maxlen
    integer :: nAppend
    integer :: rem_app
    integer :: curval, iVal, iOld, iAdd
    integer, allocatable :: newSorted(:)
    !------------------------------------------------------------------------

    if (size(val) == 0) return

    veclen = size(val)
    maxlen = veclen + me%nVals

    allocate(newSorted(maxlen))

    addedval = .false.

    iOld = 1
    iAdd = 1

    nAppend = 0
    curVal = 0

    ! Select the first entry before the loop unconditionally without checks
    ! for uniqueness (nothing to check against yet).
    if ( me%val(me%sorted(iOld)) <= val(iAdd) ) then
      curval = curval + 1
      newSorted(curval) = me%sorted(iOld)
      lastval = me%val(me%sorted(iOld))
      iOld = iOld + 1
    else
      curval = curval + 1
      nAppend = nAppend + 1
      newSorted(curval) = me%nVals + nAppend
      lastval = val(iAdd)
      if (present(pos))  pos(iAdd) = newSorted(curval)
      addedVal(iAdd) = .true.
      iAdd = iAdd + 1
    end if

    do iVal=2,maxlen

      if ( (iAdd <= veclen) .and. (iOld <= me%nVals) ) then

        if ( me%val(me%sorted(iOld)) <= val(iAdd) ) then

          ! The original list's values are appended to newsorted before
          ! the additional list is appended.
          curval = curval + 1
          newSorted(curval) = me%sorted(iOld)
          lastval = me%val(me%sorted(iOld))
          iOld = iOld + 1

        else

          ! Only append the value to unique lists, if it is not yet in the list.
          ! (If it is already in the list, it has to be the previous (curval-1)
          !  entry.)
          if ( lastval < val(iAdd) ) then
            nAppend = nAppend + 1
            curval = curval + 1
            newSorted(curval) = me%nVals + nAppend
            lastval = val(iAdd)
            addedVal(iAdd) = .true.
          end if
          if (present(pos)) pos(iAdd) = newSorted(curval)
          iAdd = iAdd + 1

        end if

      else

        ! Reached the end of one or both of the sorted lists.
        EXIT

      end if

    end do

    if (iOld <= me%nVals) then
      ! Still some values from the original list left.
      newSorted(curVal+1:me%nvals+nAppend) = me%sorted(iOld:me%nVals)
    end if

    if (iAdd <= veclen) then
      ! Still some values from the list to append left.
      rem_app = iAdd
      do i = rem_app,veclen
        if ( lastval < val(iAdd) ) then
          nAppend = nAppend + 1
          curval = curval + 1
          newSorted(curVal) = me%nVals + nAppend
          lastval = val(iAdd)
          addedVal(iAdd) = .true.
        end if
        if (present(pos)) pos(iAdd) = newSorted(curval)
        iAdd = iAdd + 1
      end do
    end if

    if (me%nVals > huge(me%nVals)-nAppend) then
       write(*,*) "Reached end of integer range for dynamic array!"
       write(*,*) "Aborting!!"
       STOP
    end if

    if (me%nVals + nAppend > me%containersize) then
      call expand( me        = me,      &
        &          increment = nAppend, &
        &          length    = length   )
    end if
    me%sorted(:me%nVals+nAppend) = newSorted(:me%nVals+nAppend)
    curval = me%nVals
    do iAdd=1,veclen
      if (addedVal(iAdd)) then
        curval = curval + 1
        me%val(curval) = val(iAdd)
      end if
    end do
    me%nVals = me%nVals + nAppend

    if( present( wasAdded ) ) wasAdded = addedval

  end subroutine append_DA_vec?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !!
  subroutine truncate_DA_?tname?(me)
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: me !< Array to sortTruncate
    !------------------------------------------------------------------------
    ?tstring?, allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------

    if (me%nVals < me%ContainerSize) then
      allocate(swpval(me%nVals))
      allocate(swpsort(me%nVals))

      swpval = me%val(:me%nVals)
      swpsort = me%sorted(:me%nVals)

      call move_alloc(swpval, me%val)
      call move_alloc(swpsort, me%sorted)

      me%containerSize = me%nVals
    end if

  end subroutine truncate_DA_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Empty all contents of the array without changing the size or status of any
  !! array
  !!
  subroutine empty_DA_?tname?(me)
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: me !< Array to sortTruncate
    !------------------------------------------------------------------------
    ! Reset the number of entries
    me%nVals = 0
  end subroutine empty_DA_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Fixing the dynamic array
  !!
  !! Truncate the array after the last valid entry and hence cut off the empty
  !! trailing empty entries
  !! Store the array in the sorted order according to the sorted( ) array
  !!
  subroutine sortTruncate_DA_?tname?(me)
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: me !< Array to sortTruncate
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: tArray !< temporary Array
    integer :: iVal
    integer :: dPos
    !------------------------------------------------------------------------
    ! Allocate the temporary array
    call init( me = tArray, length = me%nVals )
    ! Copy the entries in a sorted fashion into the temporary array
    do iVal = 1, me%nVals
      call append( me = tArray, val = me%val( me%sorted( iVal )), &
           &       pos = dPos)
    enddo
    call destroy( me = me )

    me = tArray
    call destroy( me = tArray )

  end subroutine sortTruncate_DA_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Expanding the dynamic array
  !!
  !! This is a helping subroutine, which doubles the container
  !! of the given dynamic array. As the container might be
  !! initially 0-sized, a module variable minLength has been introduced, which
  !! is used here, to at least create a container of this size.
  subroutine expand_DA_?tname?(me, increment, length)
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type) :: me !< Array to resize
    integer, optional :: increment !< Used for vector append
    !> optional length to expand the array
    integer, intent(in), optional :: length
    !------------------------------------------------------------------------
    ?tstring?, allocatable :: swpval(:)
    integer, allocatable :: swpsort(:)
    !------------------------------------------------------------------------
    integer :: addvals, expLen
    !------------------------------------------------------------------------

    addvals = 1
    if (present(increment)) addvals = increment

    if (addvals > 0) then

      ! If length is present, use that, otherwise double the size
      if( present( length ) ) then
        expLen = length
      else
        ! Set the global minimum length, if doubling would be smaller than that
        expLen = max(me%ContainerSize, minLength)
      end if

      ! Check whether all elements will fit
      if( addvals > expLen ) then
        expLen = addVals
      end if

      ! Check whether the new size will exceed the max container size.
      if( (huge(me%ContainerSize) - expLen) <= me%ContainerSize ) then
        ! if so, expand to the maximum size
        me%ContainerSize = huge(me%ContainerSize)
      else
        ! if not, expand to the calculated size
        me%ContainerSize = me%ContainerSize + expLen
      end if

      ! Only need to do something, if there are actually values to append.
      if (me%nVals > 0) then

        allocate(swpval(me%containerSize))
        swpval(1:me%nVals) = me%Val(1:me%nVals)
        call move_alloc( swpval, me%Val )

        allocate(swpsort(me%containerSize))
        swpsort(1:me%nVals) = me%sorted(1:me%nVals)
        call move_alloc( swpsort, me%sorted )

      else ! me%nVals == 0

        if( allocated(me%Val) ) &
          deallocate(me%Val)
        allocate(me%Val(me%ContainerSize))
        if( allocated(me%sorted) ) &
          deallocate(me%sorted)
        allocate(me%sorted(me%ContainerSize))

      end if

    end if

  end subroutine expand_DA_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> Return the sorted position of a value in the given dynamic array
  !!
  !! If the value was not found,
  !!  - return 0 if nextIfNotFound = .false.
  !!  - return position at the end if nextIfNotFound = .true.
  function SortPosOfVal_?tname?(me, val, nextIfNotFound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type), intent(in) :: me !< Dynamic Array
    ?argstring?, intent(in) :: val !< Value to look for
    !> Flag to indicate, if the next entry in the list should be returned,
    !! if the searched one is not found.
    logical, intent(in), optional :: nextIfNotFound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< Position of val in the sorted list, 0 if not found
    !------------------------------------------------------------------------
    logical :: retNext
    integer :: lb, ub
    integer :: mid
    ?tstring? :: lb_val, ub_val
    ?tstring? :: mid_val
    !------------------------------------------------------------------------

    retNext = .false.
    if (present(nextIfNotFound)) retNext = nextIfNotFound

    lb = 1
    ub = me%nVals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    pos = 0
    if (retNext) pos = lb

    !> Binary search on sorted list
    do while(ub >= lb)
      lb_val = me%Val(me%sorted(lb))

      if (val < lb_val) then
        if (retNext) pos = lb
        exit
      end if

      ub_val = me%Val(me%sorted(ub))

      if (val > ub_val) then
        if (retNext) pos = ub+1
        exit
      end if

      ! Safe guard against integer limit overflow
      mid = lb + (ub-lb) / 2
      mid_val = me%Val(me%sorted(mid))
      if (val == mid_val) then
        pos = mid
        exit
      end if
      if (val > mid_val) then
        lb = mid + 1
      else
        ub = mid - 1
      end if
    end do
  end function SortPosOfVal_?tname?
! ***************************************************************************** !


! ***************************************************************************** !
  !> The actual position of a given value in the dynamic array
  !!
  !! Most likely this is what you need in codes, using this
  !! data structure, it first does the binary search on the sorted
  !! values with SortPosOfVal_?tname? and then returns the looked
  !! up position in the original unsorted array, which corresponds
  !! to the position returned by the append routine.
  function PosOfVal_?tname?(me, val, nextIfNotFound, lower, upper) result(pos)
    !------------------------------------------------------------------------
    type(dyn_?tname?Array_type), intent(in) :: me !< Dynamic Array
    ?argstring?, intent(in) :: val !< Value to search for
    !> Flag to indicate, if the position of the next entry in the sorted
    !! list should be returned instead, if val is not found.
    logical, intent(in), optional :: nextIfNotFound
    integer, intent(in), optional :: lower !< lower search limit
    integer, intent(in), optional :: upper !< upper search limit
    integer :: pos !< Position in the array of the searche value, 0 if not found
    !------------------------------------------------------------------------
    integer :: sortpos
    integer :: lb, ub
    !------------------------------------------------------------------------

    lb = 1
    ub = me%nVals

    if( present( lower ) ) lb = lower
    if( present( upper ) ) ub = upper

    sortpos = SortedPosOfVal(me, val, nextIfNotFound, lb, ub)

    ! If result (sorted pos)
    if ((sortpos <= me%nVals) .and. (sortpos > 0)) then
      pos = me%sorted(sortpos)
    else
      pos = sortpos
    end if
  end function PosOfVal_?tname?
! ***************************************************************************** !
?? end text DA_impltxt
