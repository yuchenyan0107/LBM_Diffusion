! Copyright (c) 2014, 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
! Copyright (c) 2016 Tobias Schneider <tobias1.schneider@student.uni-siegen.de>
! Copyright (c) 2016, 2019 Harald Klimach <harald.klimach@uni-siegen.de>
! Copyright (c) 2016-2018 Peter Vitt <peter.vitt2@uni-siegen.de>
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
!> This UTEST is to test variable system.
!! load variables and access state and derive quantities from
!! c pointer methoddata
program tem_spacetime_fun_test
  use, intrinsic :: iso_c_binding, only: C_NEW_LINE, c_ptr, c_loc, c_f_pointer

  use env_module,               only: rk, solSpecLen, fin_env
  use tem_bc_prop_module,       only: tem_bc_prop_type
  use tem_general_module,       only: tem_general_type
  use treelmesh_module,         only: treelmesh_type
  use tem_utestEnv_module,      only: load_env
  use tem_tools_module,         only: upper_to_lower
  use tem_spacetime_fun_module, only: tem_spacetime_fun_type,             &
    &                                 tem_load_spacetime,                 &
    &                                 tem_spacetime_for
  use tem_logging_module,       only: tem_logging_init, logUnit

  use aotus_module,             only: open_config_chunk, close_config, &
    &                                 flu_state
  use aot_table_module,         only: aot_table_open, aot_table_close, &
    &                                 aot_table_length, aot_get_val

  !mpi!nprocs = 1

  implicit none

  character, parameter :: nl = C_NEW_LINE

  character(len=solSpecLen), parameter :: stfun_Conf =               &
    &    'function luafun_vector(x,y,z,t)' // nl                     &
    & // '  return {x,y,z}' // nl                                    &
    & // 'end' // nl                                                 &
    & // 'function luafun_scalar(x,y,z,t)' // nl                     &
    & // '  return x' // nl                                          &
    & // 'end' // nl                                                 &
    & // 'shape = {' // nl                                           &
    & // '  kind = "canoND",' // nl                                  &
    & // '  object = {' // nl                                        &
    & // '    origin = { 0.0, 0.5, 0.5 },' //nl                      &
    & // '    vec = { 1.0, 0.0, 0.0 }' // nl                         &
    & // '  }' // nl                                                 &
    & // '}' // nl                                                   &
    & // 'e1 = {' // nl                                              &
    & // '  kind = "canoND",' // nl                                  &
    & // '  object = { origin = { 0.25, 0.25, 0.25 } }' // nl        &
    & // '}' // nl                                                   &
    & // 'e2 = {' // nl                                              &
    & // '  kind = "canoND",' // nl                                  &
    & // '  object = { origin = { 0.75, 0.25, 0.25 } }' // nl        &
    & // '}' // nl                                                   &
    & // 'scalar_fun = luafun_scalar' // nl                          &
    & // 'scalar_const = 1.0' // nl                                  &
    & // 'vector_fun = luafun_vector' // nl                          &
    & // 'vector_const = { 1.0, 2.0, 3.0 }' // nl                    &
    & // 'vector_const_multitab = {' // nl                           &
    & // '  { 1.0, 2.0, 3.0},' // nl                                 &
    & // '  { 4.0, 5.0, 6.0}' // nl                                  &
    & // '}' // nl                                                   &
    & // 'scalar_fun_shapeall = { fun = luafun_scalar }' // nl       &
    & // 'scalar_const_shapeall = { const = scalar_const }' // nl    &
    & // 'predefined_shapeall = {' // nl                             &
    & // '  predefined = "combined",' // nl                          &
    & // '  temporal = 1.0,' // nl                                   &
    & // '  spatial = 1.0' // nl                                     &
    & // '}' // nl                                                   &
    & // 'scalar_fun_shape = {' // nl                                &
    & // '  fun = luafun_scalar,' // nl                              &
    & // '  shape = shape' // nl                                     &
    & // '}' // nl                                                   &
    & // 'scalar_const_shape = { const = 1.0, shape = shape }' // nl &
    & // 'vector_fun_shape = {' // nl                                &
    & // '  fun = luafun_vector,' // nl                              &
    & // '  shape = shape' // nl                                     &
    & // '}' // nl                                                   &
    & // 'vector_const_shape = {' // nl                              &
    & // '  const = { 1.0, 2.0, 3.0 },' // nl                        &
    & // '  shape = shape' // nl                                     &
    & // '}' // nl                                                   &
    & // 'predefined_shape = {' // nl                                &
    & // '  predefined = "combined",' // nl                          &
    & // '  temporal = 1.0,' // nl                                   &
    & // '  spatial = scalar_fun,' // nl                             &
    & // '  shape = shape' // nl                                     &
    & // '}' // nl                                                   &
    & // 'predefined_vec_shape = {' // nl                            &
    & // '  predefined = "combined",' // nl                          &
    & // '  temporal = 1.0,' // nl                                   &
    & // '  spatial = {1.0,2.0,3.0},' // nl                          &
    & // '  shape = shape' // nl                                     &
    & // '}' // nl                                                   &
    & // 'predefined_subtable_vec_shape = {' // nl                   &
    & // '  predefined = {' // nl                                    &
    & // '    "combined",' // nl                                     &
    & // '    temporal = 1.0,' // nl                                 &
    & // '    spatial = {1.0,2.0,3.0},' // nl                        &
    & // '    shape = shape' // nl                                   &
    & // '  }' // nl                                                 &
    & // '}' // nl                                                   &
    & // 'st_fun_singleTab = { scalar_fun }' // nl                   &
    & // 'st_fun_singleTab_shape = {' // nl                          &
    & // '  fun=scalar_fun,' // nl                                   &
    & // '  shape=shape' // nl                                       &
    & // '}' // nl                                                   &
    & // 'st_fun_multiTab = { { fun=scalar_fun } }' // nl            &
    & // 'st_fun_const = {' // nl                                    &
    & // '  scalar_const,' // nl                                     &
    & // '  scalar_const*2,' // nl                                   &
    & // '  scalar_const*3' // nl                                    &
    & // '}' // nl                                                   &
    & // 'st_fun_vec_const = {' // nl                                &
    & // '  vector_const,' // nl                                     &
    & // '  vector_const,' // nl                                     &
    & // '  vector_const' // nl                                      &
    & // '}' // nl                                                   &
    & // 'st_fun_2_s = {' // nl                                      &
    & // '  scalar_fun,' // nl                                       &
    & // '  scalar_fun,' // nl                                       &
    & // '  predefined_shape' // nl                                  &
    & // '}'

  type solver_type
    integer :: nComps
    integer :: nDofs
    real(kind=rk), allocatable :: state(:)
    type(treelmesh_type) :: tree
    type(tem_bc_prop_type) :: boundary
    type(tem_general_type) :: general
  end type solver_type

  type(flu_state) :: conf
  type(solver_type), target :: solver
  integer :: errorCode, iTest, i
  logical :: res(23)
  type(tem_spacetime_fun_type) :: density
  type(tem_spacetime_fun_type) :: velocity
  character(len=8) :: fi, fl

  write(*,*) 'Hello from tem_varSys_test'
  ! load utest mesh
  call load_env( tree     = solver%tree,     &
    &            boundary = solver%boundary, &
    &            general  = solver%general   )

  call tem_logging_init( level = 10, rank = 0 )

  write(*,*) 'nElems ', solver%tree%nElems
  errorCode = 0

  call open_config_chunk(L=conf, chunk=trim(stfun_conf))
  iTest = 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,     &
    &                      conf    = conf,        &
    &                      errCode = errorCode,   &
    &                      key     = 'scalar_fun' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,       &
    &                      conf    = conf,          &
    &                      errCode = errorCode,     &
    &                      key     = 'scalar_const' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,    &
    &                      conf    = conf,        &
    &                      nComp   = 3,           &
    &                      errCode = errorCode,   &
    &                      key     = 'vector_fun' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,      &
    &                      conf    = conf,          &
    &                      nComp   = 3,             &
    &                      errCode = errorCode,     &
    &                      key     = 'vector_const' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,               &
    &                      conf    = conf,                   &
    &                      nComp   = 3,                      &
    &                      errCode = errorCode,              &
    &                      key     = 'vector_const_multitab' )
  res(iTest) = errorCode == -1

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,              &
    &                      conf    = conf,                 &
    &                      errCode = errorCode,            &
    &                      key     = 'scalar_fun_shapeall' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,                &
    &                      conf    = conf,                   &
    &                      errCode = errorCode,              &
    &                      key     = 'scalar_const_shapeall' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,              &
    &                      conf    = conf,                 &
    &                      errCode = errorCode,            &
    &                      key     = 'predefined_shapeall' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,           &
    &                      conf    = conf,              &
    &                      errCode = errorCode,         &
    &                      key     = 'scalar_fun_shape' )
  res(iTest) = errorCode == 0
  write(*,*) 'shapekind ', density%geom(:)%kind

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,             &
    &                      conf    = conf,                &
    &                      errCode = errorCode,           &
    &                      key     = 'scalar_const_shape' )
  res(iTest) = errorCode == 0
  write(*,*) 'shapekind ', density%geom(:)%kind

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,          &
    &                      conf    = conf,              &
    &                      nComp   = 3,                 &
    &                      errCode = errorCode,         &
    &                      key     = 'vector_fun_shape' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,            &
    &                      conf    = conf,                &
    &                      errCode = errorCode,           &
    &                      nComp   = 3,                   &
    &                      key     = 'vector_const_shape' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = density,           &
    &                      conf    = conf,              &
    &                      errCode = errorCode,         &
    &                      key     = 'predefined_shape' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,              &
    &                      conf    = conf,                  &
    &                      ncomp   = 3,                     &
    &                      errCode = errorCode,             &
    &                      key     = 'predefined_vec_shape' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,                       &
    &                      conf    = conf,                           &
    &                      ncomp   = 3,                              &
    &                      errCode = errorCode,                      &
    &                      key     = 'predefined_subtable_vec_shape' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,    &
    &                      conf    = conf,        &
    &                      errCode = errorCode,   &
    &                      key     = 'scalar_fun' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,      &
    &                      conf    = conf,          &
    &                      nComp   = 3,             &
    &                      errCode = errorCode,     &
    &                      key     = 'vector_const' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,          &
    &                      conf    = conf,              &
    &                      errCode = errorCode,         &
    &                      key     = 'st_fun_singleTab' )
  ! This test should fail as the scalar value in a table is neither a vector
  ! nor is there one of the keywords, in this case 'const'.
  ! See scalar_const_shapeall
  res(iTest) = errorCode == -1

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,                &
    &                      conf    = conf,                    &
    &                      errCode = errorCode,               &
    &                      key     = 'st_fun_singleTab_shape' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,         &
    &                      conf    = conf,             &
    &                      errCode = errorCode,        &
    &                      key     = 'st_fun_multiTab' )
  ! A table within a table is not allowed. Therefore this one should fail
  res(iTest) = errorCode == -1

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,      &
    &                      conf    = conf,          &
    &                      errCode = errorCode,     &
    &                      key     = 'st_fun_const' )
  res(iTest) = errorCode == 0

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,          &
    &                      conf    = conf,              &
    &                      nComp   = 3,                 &
    &                      errCode = errorCode,         &
    &                      key     = 'st_fun_vec_const' )
  ! Again, tables within tables are not allowed. This time the tables are
  ! vectors, however, they need to be nested within a multiples-stfun
  res(iTest) = errorCode == -1

  iTest = iTest + 1
  write(*,"(A)") "---------------"
  write(*,"(A,I2)") "Running test ", iTest
  call tem_load_spacetime( me      = velocity,    &
    &                      conf    = conf,        &
    &                      errCode = errorCode,   &
    &                      key     = 'st_fun_2_s' )
  ! Two scalars within a table is not allowed. It is either a vector, then
  ! the notation is not corrent, or they qre two stfuns, then they need to be
  ! encapsulated by a multiples-stfun
  write(*,*) errorCode
  res(iTest) = errorCode == -1

  if ( .not. all(res) ) then
    write(logUnit(1),*) 'Error: Failed loading spacetime function calls.'
    write(fi,'(A,I2,A)') "(A,",iTest,"I3)"
    write(fl,'(A,I2,A)') "(A,",iTest,"L3)"
    write(logUnit(1),fi) 'Test no.', (i,i=1,iTest)
    write(logUnit(1),fl) 'Result  ', res
  else
    write(*,*) 'PASSED'
  end if
  call close_config(L=conf)
  call fin_env()

end program tem_spacetime_fun_test

