! Copyright (c) 2015 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
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
program tem_dynArray_test
  use env_module,           only: long_k, rk
  use tem_dyn_array_module, only: dyn_longArray_type,   &
    &                             dyn_intarray_type,    &
    &                             dyn_realarray_type,   &
    &                             dyn_labelArray_type,  &
    &                             init, append, destroy

  implicit none
  
  call check_dyn_array( valname='int_values'   )
  call check_dyn_array( valname='real_values'  )
  call check_dyn_array( valname='label_values' )
  call check_dyn_array( valname='long_values'  )

contains
  
  subroutine check_dyn_array( valname )
    ! -------------------------------------------------------------------------!
    character(len=*), intent(in) :: valname
    ! -------------------------------------------------------------------------!
    type(dyn_intarray_type)   :: da_int
    type(dyn_realarray_type)  :: da_real
    type(dyn_labelarray_type) :: da_label
    type(dyn_longArray_type)  :: da_long
    integer :: pos(5)
    integer :: iVal
    logical :: is_sorted
    ! -------------------------------------------------------------------------!
    select case (trim(valname))

    case ('int_values')
      ! initialize dynamic array of size 2
      call init( me     = da_int, &
        &        length = 2       )
      if (da_int%containersize /= 2) then
        write(*,*) 'Unexpected containersize!'
        write(*,*) 'FAILED'
      else if (da_int%nvals /= 0) then
        write(*,*) 'Unexpected number of values in list!'
        write(*,*) 'FAILED'
      end if
      ! Creating a unique list
      call append( me  = da_int, &
        &          val = 42,     &
        &          pos = pos(1)  )

      call append( me  = da_int, &
        &          val = 4711,   &
        &          pos = pos(1)  )

      call append( me  = da_int, &
        &          val = 24,     &
        &          pos = pos(1)  )

      call append( me  = da_int, &
        &          val = 32,     &
        &          pos = pos(1)  )

      call append( me  = da_int, &
        &          val = 12,     &
        &          pos = pos(1)  )

      call append( me  = da_int, &
        &          val = 2015,   &
        &          pos = pos(1)  )

      call append( me  = da_int,             &
        &          val = [1, 2, 12, 27, 27], &
        &          pos = pos                 )
      write(*,*) 'Sorted list:', da_int%val(da_int%sorted(:da_int%nVals))

      if (da_int%nVals /= 9) then

        write(*,*) 'Unexpected number of values in sorted list!'
        write(*,*) 'FAILED'

      else

        is_sorted = .true.
        do iVal=2,da_int%nVals
          if ( da_int%val(da_int%sorted(iVal-1)) >= da_int%val(da_int%sorted(iVal)) ) then
            is_sorted = .false.
            EXIT
          end if
        end do

        if (is_sorted) then
          write(*,*) 'PASSED'
        else
          write(*,*) 'Entries are not sorted!'
          write(*,*) 'FAILED'
        end if

      end if

      write(*,*) 'Destroying dynamic array now.'
      call destroy( me = da_int )
      
      if (da_int%nVals == 0 .and. da_int%containersize == 0) then
        write(*,*) 'PASSED'
      else
        write(*,*) 'Dynamic array is not empty!'
        write(*,*) 'FAILED'
      end if

    case ('real_values')
      ! initialize dynamic array of size 2
      call init( me     = da_real,  &
        &        length = 2         )
      if (da_real%containersize /= 2) then
        write(*,*) 'Unexpected containersize!'
        write(*,*) 'FAILED'
      else if (da_real%nvals /= 0) then
        write(*,*) 'Unexpected number of values in list!'
        write(*,*) 'FAILED'
      end if
      ! Creating a unique list
      call append( me  = da_real,  &
        &          val = 42._rk,   &
        &          pos = pos(1)    )

      call append( me  = da_real,  &
        &          val = 4711._rk, &
        &          pos = pos(1)    )

      call append( me  = da_real,  &
        &          val = 24._rk,   &
        &          pos = pos(1)    )

      call append( me  = da_real,  &
        &          val = 32._rk,   &
        &          pos = pos(1)    )

      call append( me  = da_real,  &
        &          val = 12._rk,   &
        &          pos = pos(1)    )

      call append( me  = da_real,  &
        &          val = 2015._rk, &
        &          pos = pos(1)    )

      call append( me  = da_real,                                &
        &          val = [1._rk, 2._rk, 12._rk, 27._rk, 27._rk], &
        &          pos = pos                                     )

      write(*,*) 'Sorted list:', da_real%val(da_real%sorted(:da_real%nVals))

      if (da_real%nVals /= 9) then

        write(*,*) 'Unexpected number of values in sorted list!'
        write(*,*) 'FAILED'

      else

        is_sorted = .true.
        do iVal=2,da_real%nVals
          if ( da_real%val(da_real%sorted(iVal-1)) >= da_real%val(da_real%sorted(iVal)) ) then
            is_sorted = .false.
            EXIT
          end if
        end do

        if (is_sorted) then
          write(*,*) 'PASSED'
        else
          write(*,*) 'Entries are not sorted!'
          write(*,*) 'FAILED'
        end if

      end if

      write(*,*) 'Destroying dynamic array now.'
      call destroy( me = da_real )
      
      if (da_real%nVals == 0 .and. da_real%containersize == 0) then
        write(*,*) 'PASSED'
      else
        write(*,*) 'Dynamic array is not empty!'
        write(*,*) 'FAILED'
      end if


    case ('label_values')
      ! initialize dynamic array of size 2
      call init( me     = da_label, &
        &        length = 2         )
      if (da_label%containersize /= 2) then
        write(*,*) 'Unexpected containersize!'
        write(*,*) 'FAILED'
      else if (da_label%nvals /= 0) then
        write(*,*) 'Unexpected number of values in list!'
        write(*,*) 'FAILED'
      end if
      ! Creating a unique list
      call append( me  = da_label,  &
        &          val = 'f',       &
        &          pos = pos(1)     )

      call append( me  = da_label,  &
        &          val = 'h',       &
        &          pos = pos(1)     )

      call append( me  = da_label,  &
        &          val = 'b',       &
        &          pos = pos(1)     )

      call append( me  = da_label,  &
        &          val = 'd',       &
        &          pos = pos(1)     )

      call append( me  = da_label,  &
        &          val = 'c',       &
        &          pos = pos(1)     )

      call append( me  = da_label,  &
        &          val = 'g',       &
        &          pos = pos(1)     )

      call append( me  = da_label,                  &
        &          val = ['a', 'b', 'c', 'e', 'e'], &
        &          pos = pos                        )
      write(*,*) 'Sorted list:', da_label%val(da_label%sorted(:da_label%nVals))

      if (da_label%nVals /= 9) then

        write(*,*) 'Unexpected number of values in sorted list!'
        write(*,*) 'FAILED'

      else

        is_sorted = .true.
        do iVal=2,da_label%nVals
          if ( da_label%val(da_label%sorted(iVal-1)) >= da_label%val(da_label%sorted(iVal)) ) then
            is_sorted = .false.
            EXIT
          end if
        end do

        if (is_sorted) then
          write(*,*) 'PASSED'
        else
          write(*,*) 'Entries are not sorted!'
          write(*,*) 'FAILED'
        end if

      end if

      write(*,*) 'Destroying dynamic array now.'
      call destroy( me = da_label )
      
      if (da_label%nVals == 0 .and. da_label%containersize == 0) then
        write(*,*) 'PASSED'
      else
        write(*,*) 'Dynamic array is not empty!'
        write(*,*) 'FAILED'
      end if


    case ('long_values')
      ! initialize dynamic array of size 2
      call init( me     = da_long,  &
        &        length = 2         )
      if (da_long%containersize /= 2) then
        write(*,*) 'Unexpected containersize!'
        write(*,*) 'FAILED'
      else if (da_long%nvals /= 0) then
        write(*,*) 'Unexpected number of values in list!'
        write(*,*) 'FAILED'
      end if
      ! Creating a unique list
      call append( me  = da_long,     &
        &          val = 42_long_k,   &
        &          pos = pos(1)       )

      call append( me  = da_long,     &
        &          val = 4711_long_k, &
        &          pos = pos(1)       )

      call append( me  = da_long,     &
        &          val = 24_long_k,   &
        &          pos = pos(1)       )

      call append( me  = da_long,     &
        &          val = 32_long_k,   &
        &          pos = pos(1)       )

      call append( me  = da_long,     &
        &          val = 12_long_k,   &
        &          pos = pos(1)       )

      call append( me  = da_long,     &
        &          val = 2015_long_k, &
        &          pos = pos(1)       )

      call append( me  = da_long,                                               &
        &          val = [1_long_k, 2_long_k, 12_long_k, 27_long_k, 27_long_k], &
        &          pos = pos                                                    )
      write(*,*) 'Sorted list:', da_long%val(da_long%sorted(:da_long%nVals))


      if (da_long%nVals /= 9) then

        write(*,*) 'Unexpected number of values in sorted list!'
        write(*,*) 'FAILED'

      else

        is_sorted = .true.
        do iVal=2,da_long%nVals
          if ( da_long%val(da_long%sorted(iVal-1)) >= da_long%val(da_long%sorted(iVal)) ) then
            is_sorted = .false.
            EXIT
          end if
        end do

        if (is_sorted) then
          write(*,*) 'PASSED'
        else
          write(*,*) 'Entries are not sorted!'
          write(*,*) 'FAILED'
        end if

      end if

      write(*,*) 'Destroying dynamic array now.'
      call destroy( me = da_long )
      
      if (da_long%nVals == 0 .and. da_long%containersize == 0) then
        write(*,*) 'PASSED'
      else
        write(*,*) 'Dynamic array is not empty!'
        write(*,*) 'FAILED'
      end if
    end select

    end subroutine check_dyn_array

end program tem_dynArray_test
