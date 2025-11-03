! Copyright (c) 2012-2013 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2012 Manuel Hasert <m.hasert@grs-sim.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2020 Neda Ebrahimi Pour <neda.epour@uni-siegen.de>
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
program tem_construction_test
  use env_module, only: long_k
  use tem_construction_module, only: tem_find_depProc_globSearch !, tem_find_BCs_fromCoarser
  use tem_element_module, only: tem_element_type
  use tem_stencil_module, only: tem_stencilHeader_type, &
    &                           d3q19_cxDir, init, append
  use tem_topology_module, only: tem_path_type, tem_pathof

  implicit none

  logical :: passed

  passed = .true.
  write(*,*) 'Running tem_construction_test...'
  call check_find_depProc
  call check_find_BCs_fromCoarser

  if (passed) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

  subroutine check_find_depProc()
    integer(kind=long_k) :: part_first(16)
    integer(kind=long_k) :: part_last(16)
    type(tem_path_type) :: pathFirst(16)
    type(tem_path_type) :: pathLast(16)
    type(tem_path_type) :: neighpath
    integer :: depProc
    integer :: nDepProcs
    integer :: iTID
    logical :: failure

    part_first=[ 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69]
    part_last =[12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]

    pathFirst     =  tem_PathOf(Part_First)
    pathLast      =  tem_PathOf(Part_Last)

    write(*,*) ' * checking tem_find_depProc...'

    write(*,*) '   checking all 64 elements on level 2, in 16 partitions:'
    failure = .false.
    do iTID=9,72
      neighpath = tem_PathOf(int(iTID, kind=long_k))
      call tem_find_depProc_globSearch(depProc, nDepProcs, &
        &                    neighPath, 1, 16, pathFirst, pathLast)
      if (nDepProcs /= 1) then
        write(*,*) '   - ERROR for treeID ', iTID, ': nDepProcs /= 1!'
        write(*,*) '     depProc: ', depProc
        failure = .true.
      else
        if (depProc /= (iTID-5)/4) then
          write(*,*) '   - ERROR for treeID ', iTID, ': wrong partition'
          write(*,*) '     is: ', depProc, '; should be:', (iTID-5)/4
          failure = .true.
        end if
      end if
    end do

    if (failure) then
      write(*,*) '   > FAILURE in tem_find_depProc single level'
      passed = .false.
    else
      write(*,*) '   > SUCCESS in tem_find_depProc single level'
    end if

    write(*,*) '   checking all 8 elements on level 1, in 16 partitions:'
    failure = .false.
    do iTID=1,8
      neighpath = tem_PathOf(int(iTID, kind=long_k))
      call tem_find_depProc_globSearch(depProc, nDepProcs, &
        &                    neighPath, 1, 16, pathFirst, pathLast)
      if (nDepProcs /= 2) then
        write(*,*) '   - ERROR for treeID ', iTID, ': nDepProcs /= 2!'
        write(*,*) '     depProc: ', depProc
        failure = .true.
      else
        if (depProc /= (iTID*2)-1) then
          write(*,*) '   - ERROR for treeID ', iTID, ': wrong first partition'
          write(*,*) '     is: ', depProc, '; should be:', (iTID*2)-1
          failure = .true.
        end if
      end if
    end do

    if (failure) then
      write(*,*) '   > FAILURE in tem_find_depProc coarse level'
      passed = .false.
    else
      write(*,*) '   > SUCCESS in tem_find_depProc coarse level'
    end if

    write(*,*) ' --------------------------------------------'
    write(*,*)

  end subroutine check_find_depProc

  subroutine check_find_BCs_fromCoarser
    type(tem_stencilHeader_type) :: d3q19_stencil
    logical :: failure

    failure = .false.
    call init( me     = d3q19_stencil, &
      &        QQN    = 18,            &
      &        QQ     = 19,            &
      &        useAll = .true.,        &
      &        nDims  = 3,             &
      &        label  = 'd3q19',       &
      &        cxDir  = d3q19_cxDir    )

    write(*,*) ' * checking tem_find_BCs_fromCoarser...'


   if (failure) then
      write(*,*) '   > FAILURE in tem_find_depProc coarse level'
      passed = .false.
    else
      write(*,*) '   > SUCCESS in tem_find_depProc coarse level'
    end if

    write(*,*) ' --------------------------------------------'
    write(*,*)

  end subroutine check_find_BCs_fromCoarser

end program tem_construction_test
