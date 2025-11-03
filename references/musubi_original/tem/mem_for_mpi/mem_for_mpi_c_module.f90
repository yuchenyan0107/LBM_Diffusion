! Copyright (c) 2016 Harald Klimach <harald.klimach@uni-siegen.de>
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
module mem_for_mpi_module
  use mpi
  use, intrinsic :: iso_c_binding

  implicit none

  private

  interface
    subroutine alloc_mpi_mem(asize, baseptr, ierr) bind(c, name='alloc_mpi_mem')
      use, intrinsic :: iso_c_binding
      integer(kind=c_int64_t), value :: asize
      type(c_ptr) :: baseptr
      integer(kind=c_int) :: ierr
    end subroutine alloc_mpi_mem

    function free_mpi_mem(baseptr) bind(c, name='free_mpi_mem')
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: baseptr
      integer(kind=c_int) :: free_mpi_mem
    end function free_mpi_mem
  end interface

  public :: alloc_mpif_mem, free_mpif_mem


contains


  subroutine alloc_mpif_mem(asize, baseptr, ierr)
    integer(kind=MPI_ADDRESS_KIND), intent(in) :: asize
    type(c_ptr), intent(out) :: baseptr
    integer, optional, intent(out) :: ierr

    integer(kind=c_int) :: cerr

    call alloc_mpi_mem( asize   = int(asize, kind=c_int64_t), &
      &                 baseptr = baseptr,                    &
      &                 ierr    = cerr                        )

    if (present(ierr)) ierr = int(cerr)

  end subroutine alloc_mpif_mem


  subroutine free_mpif_mem(baseptr, ierr)
    type(c_ptr)       :: baseptr
    integer, optional :: ierr

    integer(kind=c_int) :: cerr

    cerr = free_mpi_mem(baseptr)
    if (present(ierr)) ierr = int(cerr)
  end subroutine free_mpif_mem

end module mem_for_mpi_module
