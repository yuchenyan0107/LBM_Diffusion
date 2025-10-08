! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
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
program tem_growArray_test
  use env_module,             only: long_k, rk, labelLen
  use tem_float_module,       only: operator(.fne.)
  use tem_grow_array_module,  only: grw_longarray_type,    &
    &                               grw_intarray_type,     &
    &                               grw_realarray_type,    &
    &                               grw_dtint2darray_type, &
    &                               grw_labelarray_type,   &
    &                               grw_logicalarray_type, &
    &                               grw_chararray_type,    &
    &                               intarray2d_type,       &
    &                               init, placeat, expand, &
    &                               truncate, append,      &
    &                               destroy, empty

  implicit none

  call check_grw_array( valname='2dint_values' )
  call check_grw_array( valname='long_values'  )
  call check_grw_array( valname='int_values'   )
  call check_grw_array( valname='real_values'  )
  call check_grw_array( valname='label_values' )
  call check_grw_array( valname='logic_values' )
  call check_grw_array( valname='char_values'  )

  contains

    subroutine check_grw_array( valname )
      ! -------------------------------------------------------------------------!
      character(len=*), intent(in) :: valname
      ! -------------------------------------------------------------------------!
      type(grw_dtint2darray_type) :: ga_2dint
      type(grw_longarray_type)    :: ga_long
      type(grw_intarray_type )    :: ga_int
      type(grw_realarray_type)    :: ga_real
      type(grw_labelarray_type)   :: ga_label
      type(grw_logicalarray_type) :: ga_logical
      type(grw_chararray_type)    :: ga_char
      type(intArray2d_type)       :: val2d
      character(len=labelLen),allocatable :: label(:)
      logical :: logic
      ! -------------------------------------------------------------------------!
      select case (trim(valname))

        case ('char_values')
          ! init growing array with size 2
          call init( me     = ga_char, &
            &        length = 2        )
          if (ga_char%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place character at pos 1
          call placeat( me  = ga_char, &
            &           val = 'ho', &  
            &           pos = 1)
          if (ga_char%val(1) /= 'h') then
            write(*,*) 'Unexpected character after placeat!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_char, &
            &          length = 3        )
          if (ga_char%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_char%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place ['h', 'e', 'l', 'l', 'o'] at pos (2:6)
          call placeat( me  = ga_char,                   &
            &           val = ['h', 'e', 'l', 'l', 'o'], &
            &           pos = 2                                )
          if (ga_char%val(3) /= 'e') then
            write(*,*) 'Unexpected character after placeat!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
        !  ! append growing array with 'octree'
        !  call append( me  = ga_logical, &
        !    &          val = .FALSE.     )
        !  if (ga_logical%val(6) .NEQV. .FALSE.) then
        !    write(*,*) 'Unexpected logical after append!'
        !    write(*,*) 'Expected .FALSE. at pos 6, got ', ga_logical%val(6)
        !    write(*,*) 'FAILED'
        !  else
        !    write(*,*) 'PASSED'
        !  end if
        !  ! truncate growing array
        !  call truncate( me = ga_logical )
        !  if (ga_logical%nvals /= ga_logical%containersize) then
        !    write(*,*) 'Containersize /= nvals after truncate!'
        !    write(*,*) 'FAILED'
        !  else
        !    write(*,*) 'PASSED'
        !  end if
        !  ! empty growing array
        !  call empty( me = ga_logical )
        !  if (ga_logical%nvals /= 0) then
        !    write(*,*) 'Unexpected nvals after empty!'
        !    write(*,*) 'Expected 0, got ', ga_logical%nvals
        !    write(*,*) 'FAILED'
        !  else
        !    write(*,*) 'PASSED'
        !  end if
        !  ! destroy growing array
        !  call destroy( me = ga_logical )
        !  if (ga_logical%containersize /= 0) then
        !    write(*,*) 'Unexpected containersize after destroy!'
        !    write(*,*) 'Expected 0, got ', ga_logical%containersize
        !    write(*,*) 'FAILED'
        !  elseif (ga_logical%nvals /= 0) then
        !    write(*,*) 'Unexpected nvals after destroy!'
        !    write(*,*) 'Expected 0, got ', ga_logical%nvals
        !    write(*,*) 'FAILED'
        !  else
        !    write(*,*) 'PASSED'
        !  end if

        case ('logical_values')
          ! init growing array with size 2
          call init( me     = ga_logical, &
            &        length = 2         )
          if (ga_logical%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place logical at pos 1
          logic = .TRUE.
          call placeat( me  = ga_logical, &
            &           val = logic,    &  
            &           pos = 1)
          if (ga_logical%val(1) .NEQV. .TRUE.) then
            write(*,*) 'Unexpected logical after placeat!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_logical, &
            &          length = 3         )
          if (ga_logical%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_logical%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place ['.TRUE.', '.FALSE.', '.TRUE.'] at pos (3:5)
          call placeat( me  = ga_logical,                      &
            &           val = [.TRUE., .FALSE., .TRUE.], &
            &           pos = 3                                )
          if (ga_logical%val(3) .NEQV. .TRUE.) then
            write(*,*) 'Unexpected logical after placeat!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! append growing array with 'octree'
          call append( me  = ga_logical, &
            &          val = .FALSE.     )
          if (ga_logical%val(6) .NEQV. .FALSE.) then
            write(*,*) 'Unexpected logical after append!'
            write(*,*) 'Expected .FALSE. at pos 6, got ', ga_logical%val(6)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! truncate growing array
          call truncate( me = ga_logical )
          if (ga_logical%nvals /= ga_logical%containersize) then
            write(*,*) 'Containersize /= nvals after truncate!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! empty growing array
          call empty( me = ga_logical )
          if (ga_logical%nvals /= 0) then
            write(*,*) 'Unexpected nvals after empty!'
            write(*,*) 'Expected 0, got ', ga_logical%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! destroy growing array
          call destroy( me = ga_logical )
          if (ga_logical%containersize /= 0) then
            write(*,*) 'Unexpected containersize after destroy!'
            write(*,*) 'Expected 0, got ', ga_logical%containersize
            write(*,*) 'FAILED'
          elseif (ga_logical%nvals /= 0) then
            write(*,*) 'Unexpected nvals after destroy!'
            write(*,*) 'Expected 0, got ', ga_logical%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if

        case ('label_values')
          ! init growing array with size 2
          call init( me     = ga_label, &
            &        length = 2         )
          if (ga_label%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place hellotest at pos 2
          allocate(label(1))
          label = 'hellotest'
          call placeat( me  = ga_label, &
            &           val = label,    &  
            &           pos = 1)
          deallocate(label)
          if (ga_label%val(1) /= 'hellotest') then
            write(*,*) 'Unexpected label after placeat!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_label, &
            &          length = 3         )
          if (ga_label%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_label%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place ['tr', 'ee', 'lm'] at pos (3:5)
          allocate(label(3))
          label = ['tr', 'ee', 'lm']
          call placeat( me  = ga_label, &
            &           val = label,    &
            &           pos = 3         )
          deallocate(label)
          if (ga_label%val(3) /= 'tr') then
            write(*,*) 'Unexpected label after placeat!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! append growing array with 'octree'
          allocate(label(1))
          label = 'octree'
          call append( me  = ga_label, &
            &          val = label     )
          if (ga_label%val(6) /= 'octree') then
            write(*,*) 'Unexpected label after append!'
            write(*,*) 'Expected "octree" at pos 6, got ', ga_label%val(6)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! truncate growing array
          call truncate( me = ga_label )
          if (ga_label%nvals /= ga_label%containersize) then
            write(*,*) 'Containersize /= nvals after truncate!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! empty growing array
          call empty( me = ga_label )
          if (ga_label%nvals /= 0) then
            write(*,*) 'Unexpected nvals after empty!'
            write(*,*) 'Expected 0, got ', ga_label%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! destroy growing array
          call destroy( me = ga_label )
          if (ga_label%containersize /= 0) then
            write(*,*) 'Unexpected containersize after destroy!'
            write(*,*) 'Expected 0, got ', ga_label%containersize
            write(*,*) 'FAILED'
          elseif (ga_label%nvals /= 0) then
            write(*,*) 'Unexpected nvals after destroy!'
            write(*,*) 'Expected 0, got ', ga_label%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if

        case ('2dint_values')
          ! init growing array with size 2
          call init( me     = ga_2dint, &
            &        length = 2         )
          if (ga_2dint%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place 1013 at pos 2
          allocate(val2d%val(1,1) )
          val2d%val(1,1) = 1013
          call placeat( me  = ga_2dint, &
            &           val = val2d,    &  
            &           pos = 2         )
          deallocate(val2d%val)
          if (ga_2dint%val(2)%val(1,1) /= 1013) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 1013, got ', ga_2dint%val(2)%val(1,1)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_2dint, &
            &          length = 3        )
          if (ga_2dint%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_2dint%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place [1, 2, 3] at pos 1
          allocate(val2d%val(1,3) )
          val2d%val(1,1) = 1
          val2d%val(1,2) = 2
          val2d%val(1,3) = 3
          call placeat( me  = ga_2dint, &
            &           val = val2d,    &
            &           pos = 1         )
          deallocate(val2d%val)
          if (ga_2dint%val(1)%val(1,3) /= 3) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 3, got ', ga_2dint%val(1)%val(1,3)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! append growing array with 4, ga_2dint%val is an allocatable type
          allocate(val2d%val(1,1) )
          val2d%val(1,1) = 4
          call append( me   = ga_2dint, &
            &          val  = val2d )
          deallocate(val2d%val)
          if (ga_2dint%val(3)%val(1,1) /= 4) then
            write(*,*) 'Unexpected value after append!'
            write(*,*) 'Expected 4 at pos 3, got ', ga_2dint%val(3)%val(1,1)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! truncate growing array
          call truncate( me = ga_2dint )
          if (ga_2dint%nvals /= ga_2dint%containersize) then
            write(*,*) 'Containersize /= nvals after truncate!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! empty growing array
          call empty( me = ga_2dint )
          if (ga_2dint%nvals /= 0) then
            write(*,*) 'Unexpected nvals after empty!'
            write(*,*) 'Expected 0, got ', ga_2dint%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! destroy growing array
          call destroy( me = ga_2dint )
          if (ga_2dint%containersize /= 0) then
            write(*,*) 'Unexpected containersize after destroy!'
            write(*,*) 'Expected 0, got ', ga_2dint%containersize
            write(*,*) 'FAILED'
          elseif (ga_2dint%nvals /= 0) then
            write(*,*) 'Unexpected nvals after destroy!'
            write(*,*) 'Expected 0, got ', ga_2dint%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if

        case ('long_values')
          ! init growing array with size 2
          call init( me     = ga_long, &
            &        length = 2        )
          if (ga_long%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place 1013_long_k at pos 2
          call placeat( me  = ga_long,      &
            &           val = 1013_long_k,  &
            &           pos = 2             )
          if (ga_long%val(2) /= 1013_long_k) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 1013_long_k, got ', ga_long%val(2)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_long, &
            &          length = 3        )
          if (ga_long%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_long%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place [1_long_k, 2_long_k, 3_long_k] at pos 1
          call placeat( me  = ga_long,                        &
            &           val = [1_long_k, 2_long_k, 3_long_k], &
            &           pos = 1                               )
          if (ga_long%val(3) /= 3_long_k) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 3_long_k at pos 3, got ', ga_long%val(3)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! append growing array with 4_long_k
          call append( me  = ga_long,                      &
            &          val = [4_long_k,5_long_k,6_long_k]  )
          if (ga_long%val(6) /= 6_long_k) then
            write(*,*) 'Unexpected value after append!'
            write(*,*) 'Expected 6_long_k at pos 6, got ', ga_long%val(6)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! truncate growing array
          call truncate( me = ga_long )
          if (ga_long%nvals /= ga_long%containersize) then
            write(*,*) 'Containersize /= nvals after truncate!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! empty growing array
          call empty( me = ga_long )
          if (ga_long%nvals /= 0) then
            write(*,*) 'Unexpected nvals after empty!'
            write(*,*) 'Expected 0, got ', ga_long%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! destroy growing array
          call destroy( me = ga_long )
          if (ga_long%containersize /= 0) then
            write(*,*) 'Unexpected containersize after destroy!'
            write(*,*) 'Expected 0, got ', ga_long%containersize
            write(*,*) 'FAILED'
          elseif (ga_long%nvals /= 0) then
            write(*,*) 'Unexpected nvals after destroy!'
            write(*,*) 'Expected 0, got ', ga_long%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if

        case ('int_values')
          ! init growing array with size 2
          call init( me     = ga_int, &
            &        length = 2        )
          if (ga_int%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place 1013 at pos 2
          call placeat( me  = ga_int, &
            &           val = 1013,   &
            &           pos = 2       )
          if (ga_int%val(2) /= 1013) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 1013, got ', ga_int%val(2)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_int, &
            &          length = 3       )
          if (ga_int%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_int%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place [1, 2, 3] at pos 1
          call placeat( me  = ga_int,    &
            &           val = [1, 2, 3], &
            &           pos = 1          )
          if (ga_int%val(3) /= 3) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 3 at pos 3, got ', ga_int%val(3)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! append growing array with 4
          call append( me  = ga_int, &
            &          val = [4,5,6] )
          if (ga_int%val(6) /= 6) then
            write(*,*) 'Unexpected value after append!'
            write(*,*) 'Expected 6 at pos 6, got ', ga_int%val(6)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! truncate growing array
          call truncate( me = ga_int )
          if (ga_int%nvals /= ga_int%containersize) then
            write(*,*) 'Containersize /= nvals after truncate!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! empty growing array
          call empty( me = ga_int )
          if (ga_int%nvals /= 0) then
            write(*,*) 'Unexpected nvals after empty!'
            write(*,*) 'Expected 0, got ', ga_int%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! destroy growing array
          call destroy( me = ga_int )
          if (ga_int%containersize /= 0) then
            write(*,*) 'Unexpected containersize after destroy!'
            write(*,*) 'Expected 0, got ', ga_int%containersize
            write(*,*) 'FAILED'
          elseif (ga_int%nvals /= 0) then
            write(*,*) 'Unexpected nvals after destroy!'
            write(*,*) 'Expected 0, got ', ga_int%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if

        case ('real_values')
          ! init growing array with size 2
          call init( me     = ga_real, &
            &        length = 2        )
          if (ga_real%containersize /= 2) then
            write(*,*) 'Unexpected containersize after init!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place 1013._rk at pos 2
          call placeat( me  = ga_real,  &
            &           val = 1013._rk, &
            &           pos = 2         )
          if (ga_real%val(2) .fne. 1013._rk) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 1013._rk, got ', ga_real%val(2)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! expand growing array with length 3
          call expand( me     = ga_real, &
            &          length = 3        )
          if (ga_real%containersize /= 5) then
            write(*,*) 'Unexpected containersize after expand!'
            write(*,*) 'Expected 5, got ', ga_real%containersize
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! place [1._rk, 2._rk, 3._rk] at pos 1
          call placeat( me  = ga_real,               &
            &           val = [1._rk, 2._rk, 3._rk], &
            &           pos = 1                      )
          if (ga_real%val(3) .fne. 3._rk) then
            write(*,*) 'Unexpected value after placeat!'
            write(*,*) 'Expected 3._rk at pos 3, got ', ga_int%val(3)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! append growing array with 4._rk
          call append( me  = ga_real,            &
            &          val = [4._rk,5._rk,6._rk] )
          if (ga_real%val(6) .fne. 6._rk) then
            write(*,*) 'Unexpected value after append!'
            write(*,*) 'Expected 6._rk at pos 6, got ', ga_real%val(6)
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! truncate growing array
          call truncate( me = ga_real )
          if (ga_real%nvals /= ga_real%containersize) then
            write(*,*) 'Containersize /= nvals after truncate!'
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! empty growing array
          call empty( me = ga_real )
          if (ga_real%nvals /= 0) then
            write(*,*) 'Unexpected nvals after empty!'
            write(*,*) 'Expected 0, got ', ga_real%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
          ! destroy growing array
          call destroy( me = ga_real )
          if (ga_real%containersize /= 0) then
            write(*,*) 'Unexpected containersize after destroy!'
            write(*,*) 'Expected 0, got ', ga_real%containersize
            write(*,*) 'FAILED'
          elseif (ga_real%nvals /= 0) then
            write(*,*) 'Unexpected nvals after destroy!'
            write(*,*) 'Expected 0, got ', ga_real%nvals
            write(*,*) 'FAILED'
          else
            write(*,*) 'PASSED'
          end if
      end select

    end subroutine check_grw_array
    
end program tem_growArray_test
