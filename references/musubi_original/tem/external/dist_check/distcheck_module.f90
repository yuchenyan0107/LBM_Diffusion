module distCheck_module
  use, intrinsic :: iso_c_binding

  implicit none

  private

  integer, parameter :: dk = selected_real_kind(15)

  public :: checksum

  interface checksum
    module procedure checksum_sk1, checksum_sk2, checksum_sk3, checksum_sk4, &
      &              checksum_sk5
    module procedure checksum_dk1, checksum_dk2, checksum_dk3, checksum_dk4, &
      &              checksum_dk5
  end interface

  interface
    subroutine distCRC(dat, length, comm, buf) bind(c, name="distCRC")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: dat
      integer(kind=c_long), value :: length
      integer, value :: comm
      type(c_ptr), value :: buf
    end subroutine distCRC
  end interface

contains

  function cs_string(dat, length, comm) result(cs)
    type(c_ptr) :: dat
    integer(kind=c_long), intent(in) :: length
    integer, intent(in) :: comm
    character(len=8) :: cs

    character, target :: string(9)
    integer :: i

    call distCRC(dat, length, comm, c_loc(string))
    do i=1,8
      cs(i:i) = string(i)
    end do

  end function cs_string

  function checksum_sk1(dat, length, comm) result(check)
    integer, intent(in) :: length(1)
    real(kind=c_float), intent(in), target :: dat(length(1))
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*4

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_sk1


  function checksum_sk2(dat, length, comm) result(check)
    integer, intent(in) :: length(2)
    real(kind=c_float), intent(in), target :: dat(length(1), length(2))
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*4

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_sk2


  function checksum_sk3(dat, length, comm) result(check)
    integer, intent(in) :: length(3)
    real(kind=c_float), intent(in), target :: dat( length(1), &
      &                                            length(2), &
      &                                            length(3)  )
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*length(3)*4

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_sk3


  function checksum_sk4(dat, length, comm) result(check)
    integer, intent(in) :: length(4)
    real(kind=c_float), intent(in), target :: dat( length(1), &
      &                                            length(2), &
      &                                            length(3), &
      &                                            length(4)  )
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*length(3)*length(4)*4

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_sk4


  function checksum_sk5(dat, length, comm) result(check)
    integer, intent(in) :: length(5)
    real(kind=c_float), intent(in), target :: dat( length(1), &
      &                                            length(2), &
      &                                            length(3), &
      &                                            length(4), &
      &                                            length(5)  )
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*length(3)*length(4)*length(5)*4

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_sk5


  function checksum_dk1(dat, length, comm) result(check)
    integer, intent(in) :: length(1)
    real(kind=c_double), intent(in), target :: dat(length(1))
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*8

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_dk1


  function checksum_dk2(dat, length, comm) result(check)
    integer, intent(in) :: length(2)
    real(kind=c_double), intent(in), target :: dat(length(1), length(2))
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*8

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_dk2


  function checksum_dk3(dat, length, comm) result(check)
    integer, intent(in) :: length(3)
    real(kind=c_double), intent(in), target :: dat( length(1), &
      &                                            length(2), &
      &                                            length(3)  )
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*length(3)*8

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_dk3


  function checksum_dk4(dat, length, comm) result(check)
    integer, intent(in) :: length(4)
    real(kind=c_double), intent(in), target :: dat( length(1), &
      &                                            length(2), &
      &                                            length(3), &
      &                                            length(4)  )
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*length(3)*length(4)*8

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_dk4


  function checksum_dk5(dat, length, comm) result(check)
    integer, intent(in) :: length(5)
    real(kind=c_double), intent(in), target :: dat( length(1), &
      &                                            length(2), &
      &                                            length(3), &
      &                                            length(4), &
      &                                            length(5)  )
    integer, intent(in) :: comm
    character(len=8) :: check

    integer(kind=c_long) :: nBytes

    nBytes = length(1)*length(2)*length(3)*length(4)*length(5)*8

    check = cs_string(c_loc(dat), nBytes, comm)

  end function checksum_dk5

end module distCheck_module
