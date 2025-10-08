! bof
! **********************************************************************
! Fortran 95 program coco

! **********************************************************************
! Source Control Strings

! $Id: coco.f90,v 1.30 2007/06/25 19:08:22 dan Exp dan $

! **********************************************************************
!  Copyright 2003 Purple Sage Computing Solutions, Inc.
!  All Rights Reserved

!   This program is free software; you can redistribute it and/or
!   modify it under the terms of the GNU General Public
!   License as published by the Free Software Foundation; either
!   version 2 of the License, or (at your option) any later version.

!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!   General Public License for more details.

!   You should have received a copy of the GNU General Public
!   License along with this program; if not, write to the Free
!   Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

! To report bugs, suggest enhancements, etc. to the Authors,
! Contact:
!    Purple Sage Computing Solutions, Inc.
!                               send email to dnagle@erols.com
!                                   or fax to 703 471 0684 (USA)
!                                  or mail to 12142 Purple Sage Ct.
!                                             Reston, VA 20194-5621 USA

! **********************************************************************
! coco implements Part 3: Conditional Compilation

! **********************************************************************

!  coco reads

!     input source file(s)- named via command line, or stdin
!     named set file- name taken from the output filename, or coco.set

!  coco writes

!     output source file- named via command line, or stdout
!     logfile- named via ??logfile directive, or stderr

!  coco temp files

!     scratch- hold the contents of text blocks while counting their size

!  coco uses

!     f2kcli (http://www.winteracter.com/f2kcli)

!  coco constants

!     coco_rcs_id- this file's rcs id string
!     *_unit- logical unit numbers
!     *_fmt- formats
!     *_len- lengths of character entities
!     alter_*- the alter states
!     if_*- the current state of if-directive processing

!  coco types

!     file_t- filename and unit
!     path_t- include directory
!     symbol_t- integer, logical, macro, text symbol
!     if_t- if block and state
!     state_t- set of options
!     report_t- statistics

!  coco data

!  coco library

!     add_directory() appends a directory to the path list
!     add_integer() build an integer entry and store it on the symbol list
!     add_logical() build a logical entry and store it on the symbol list
!     add_macro() build a macro entry and store it on the symbol list
!     add_symbol() append an entry to the symbol list
!     add_text() build a text entry and store it on the symbol list
!     blank_compress_lower_case() read lines return statement for processing
!     close_file() closes a named file and unlinks pointers
!     close_scratch() closes a scratch file
!     coco() main program
!     copy_set_file() appends the set file to the end of the output
!     edit_coco_strings() edit ?coco? in source
!     edit_date_time_strings() edit ?date? and ?time? strings in source lines
!     edit_file_line_strings() edit ?file? and ?line? strings in source lines
!     edit_integer_strings() edit integer substitutions
!     edit_logical_strings() edit logical substitutions
!     edit_macro_strings() edit macro substitutions
!     edit_source_line() edit a source line with substitutions
!     eval_int_expr() return an integer value from an expression
!     eval_int_primary() return a value from an integer primary
!     eval_log_expr() return a logical value from an expression
!     eval_log_primary() return a value from a logical primary
!     eval_rel_expr() return a logical value from a comparison
!     gather_coco_statement() signal when a complete statement has been read
!     get_integer_value() return the value of an integer symbol
!     get_logical_value() return the value of a logical value
!     get_macro_name() return a macro name from statement
!     get_next_integer() return a pointer to the next integer on the symbol list
!     get_next_logical() return a pointer to the next logical on the symbol list
!     get_next_macro() return a pointer to the next macro on the symbol list
!     get_symbol_name() return integer or logical name from statement
!     get_text_name() return a text block name from statement
!     get_text_ptr() return pointer to a text block name
!     getopt() get each command line word
!     initialize_coco() set derived types variables and pointers to initial state
!     integer_or_logical() determine whether an expression is integer or logical
!     is_coco_statement() decide whether line is a directive or a comment
!     msg_continue() log message and continue
!     msg_quit() log message and stop
!     open_file() opens a file and sets pointers
!     open_scratch() opens a scratch file
!     print_help() print command line options to stderr
!     process_actual_arglist() edit macro definition dummy args
!     process_alter_directive() set alter state from directive
!     process_alter_option() set alter state from -a
!     process_assert_directive() process a coco assert statement
!     process_assertif_directive() process a coco assertif statement
!     process_block_directive() process a coco statement within a text block
!     process_coco_statement() process a coco directive in a source file
!     process_command_line() process command line options including filenames
!     process_copy_directive() process a copy statement
!     process_copyif_directive() process a copyif statement
!     process_directory_directive() process a setfile directory directive
!     process_doc_directive() processes a doc directive
!     process_docfile_directive() processes a docfile directive
!     process_dummy_arglist() process a macro or text definition dummy arg list
!     process_dump_directive() process a dump directive
!     process_edit_directive() process a setfile edit directive
!     process_else_directive() process a coco else statement
!     process_elseif_directive() process a coco else if statement
!     process_endif_directive() process a coco end if statement
!     process_if_directive() process a coco if statement
!     process_ifdef_directive() true if symbol is defined
!     process_ifndef_directive() true if symbol is not defined
!     process_include_directive() process a source file include directive
!     process_include_option() set include directory from -Idir
!     process_input_file() read a source file and process it
!     process_integer_assignment() assign a value to an integer symbol
!     process_integer_constant() define an integer constant
!     process_integer_declaration() process an integer statement
!     process_logfile_directive() process a setfile logfile directive
!     process_logical_assignment() assign a value to a logical symbol
!     process_logical_constant() define a logical constant
!     process_logical_declaration() process a logical statement
!     process_macro_declaration() process a macro statement
!     process_message_directive() process a coco message statement
!     process_number_directive() set numbering from setfile directive
!     process_output_directive() opens a new output file
!     process_set_file() read and process the setfile
!     process_set_statement() process one setfile statement
!     process_stop_directive() process a coco stop statement
!     process_summary_directive() process a setfile summary directive
!     process_symbol_option() define symbol from -Dname
!     process_text_directive() process a begin text statement
!     process_undefine_directive() processes an undefine directive
!     process_verbose_directive() process a setfile verbose directive
!     process_warn_directive() process a setfile warn directive
!     process_wrap_value() process the length supplied on a -w or ??wrap
!     remove_symbol() remove a symbol entry from the symbol list
!     replace_substring() globally replace a substring in a string
!     seek_close_paren() return the location of the matching parenthesis
!     seek_directory() seek an directory to find an include file
!     seek_include_file() find an include file
!     seek_log_primary() find the next logical operator
!     seek_set_file() try to open setfile
!     seek_symbol_name() find an integer or logical name
!     set_option_defaults() check options after command line and setfile processing
!     unquote_string() return an unquoted string
!     verify_actual_args() check actual arguments for need of parenthesis
!     verify_dummy_args() verify macro or text dummy args
!     verify_macro_value() verify a macro's value
!     verfiy_text_directive() verify which directives are in a text block
!     wrap_source_line() ensure that a source line exceeds not length
!     write_coco_line() write a coco line
!     write_options() write the current options to the logfile
!     write_report() log summary statistics
!     write_source_line() write a source line

! **********************************************************************

!  coco

! **********************************************************************

program coco

!  coco implements ISO/IEC 1539-3 Conditional Compilation standard

!  coco steps

!  1. call process_command_line() to read command line, get filenames & options
!  2. call process_setfile() to read the setfile, if there is one
!  3. open the output file, if named, use stdout if not
!  4. open the input file(s), if named, use stdin if not
!  5. call process_input_file() to process the input file(s) & write the output file
!  6. copy the setfile contents to the output file
!  7. close all files
!  8. call write_report() print summary information

! **********************************************************************

!  coco uses modules

! **********************************************************************

!  Winteracter Fortran 2003 command line access module
!  http://www.winteracter.com/f2kcli

!  use f2kcli

! **********************************************************************

!  explicit declaration

implicit none

! **********************************************************************

!  coco RCS strings

! **********************************************************************

!  program source filename supplied by RCS

character( len= *), parameter :: coco_rcs_id = &
   '$Id: coco.f90,v 1.30 2007/06/25 19:08:22 dan Exp dan $'

! **********************************************************************

!  coco constants

! **********************************************************************

!  coco logical unit numbers

!  scheme for logical unit numbers:

!  The setfile is read first and processed.
!  The setfile is then closed.  The output is opened
!  (if need be), and then the input file is opened (again, if need be).
!  If an include file is encountered, the logical unit numbers used are
!  computed by adding to the current input unit number.  If the current
!  input file is stdin, read_unit is used first, then successive unit
!  numbers are used for nested include files.  When all input files have
!  been read, the set scratch file is copied to the output file.  All
!  Fortran files are then closed.  The summary is written to the output.
!  A text block is copied to the text_unit to count the number of lines.
!  Then it is copied back to a character array and the text_unit is closed.

! **********************************************************************

!  global constants

!  the unit names are provided in the Fortran 2003 intrinsic module iso_fortran_env

integer, parameter :: input_unit = 5
integer, parameter :: output_unit = 6
integer, parameter :: error_unit = 0

!  logfile unit else use error_unit, +4 tries to avoid plot_unit, punch_unit, etc.

integer, parameter :: log_unit = max( input_unit, output_unit, error_unit) + 4

!  documentation unit for doc ... end doc lines

integer, parameter :: doc_unit = log_unit + 1

!  scratch unit for text scratch files

integer, parameter :: text_unit = doc_unit + 1

!  setfile unit

integer, parameter :: set_unit = text_unit + 1

!  output unit if named output file (else use unit= *)

integer, parameter :: write_unit = set_unit + 1

!  input unit if named input file (else use unit= *)

integer, parameter :: read_unit = write_unit + 1

! **********************************************************************

!  coco formats

! **********************************************************************

!  used to read/write lines

character( len= *), parameter :: string_fmt = '( a)'

!  used to write reports

character( len= *), parameter :: integer_fmt = '( a, i10)'

!  used to write reports

character( len= *), parameter :: directory_fmt = '( a, i0, a)'

!  length of format strings

integer, parameter :: format_len = max( len( string_fmt), len( integer_fmt) )

!  length of input/output specifier strings

integer, parameter :: io_specifier_len = 16

! ----------------------------------------------------------------------

!  length of strings used to convert between integers and characters

integer, parameter :: conversion_len = 10

!  format used to convert between integers and characters

character( len= *), parameter :: conversion_fmt = '(i10)'

! **********************************************************************

!  coco character lengths

! **********************************************************************

!  these are the lengths of strings used throughout coco

! ----------------------------------------------------------------------

!  length of character storing a constant or variable name

integer, parameter :: symbol_name_len = 32

!  length of a Fortran source line

integer, parameter :: free_format_len = 132

integer, parameter :: fixed_format_len = 72

!  length used to write lines is free_format_len + len( '!?>') + len( blank)

integer, parameter :: source_line_len = free_format_len + len( '!?>') + len( ' ')

!  length of character storing filenames

integer, parameter :: filename_len = 256

!  length of character line buffers (allows for max_continuations number of continuations)

integer, parameter :: max_continuations = 39

!  buffer a whole coco statement and always have a blank at the end

integer, parameter :: buffer_len = ( max_continuations + 1) * source_line_len

! **********************************************************************

!  this string is used to initialize character variables

! ----------------------------------------------------------------------

!  null string

character( len= *), parameter :: null_string = ''

!  mark beginning of the setfile in the output

character( len= *), parameter :: mark_set_file = &
   '?? This was produced using the following SET file'

! ----------------------------------------------------------------------

!  names must be made of alphanumeric characters only

character( len= *), parameter :: alpha_chars = 'abcdefghijklmnopqrstuvwxyz'

character( len= *), parameter :: digit_chars = '0123456789'

character( len= *), parameter :: underscore = '_'

character( len= *), parameter :: alphanum_chars =  alpha_chars // digit_chars // underscore

! ----------------------------------------------------------------------

!  ascii characters change case

integer, parameter :: change_case = 32

! **********************************************************************

!  coco directives constants

! **********************************************************************

!  many character string constants' lengths are used to count past
!  the string as coco processes each statement

!  coco line and statement syntax uses the next set of character constants

! **********************************************************************

!  . separates filenames from extensions, delimits logical operators & literals

character( len= *), parameter :: dot = '.'

! ----------------------------------------------------------------------

!  constants defining coco directives, comments, separators, etc.

! ----------------------------------------------------------------------

!  coco line key ??coco_directive

character( len= *), parameter :: coco_key = '??'

!  substitution key

character( len= *), parameter :: arg_key = '?'

!  length of ?name?

integer, parameter :: target_len = len( arg_key) + symbol_name_len + len( arg_key)

!  continuation character

character( len= *), parameter :: continuation = '&'

!  blank character

character( len= *), parameter :: blank = ' '

!  tab character

character( len= *), parameter :: tab = achar( 9)

!  whitespace is blank or tab

character( len= *), parameter :: white_space = blank // tab

!  coco comment initializer

character( len= *), parameter :: comment = '!'

!  separates items within a list

character( len= *), parameter :: comma = ','

!  quotes

character( len= *), parameter :: single_quote = "'"

character( len= *), parameter :: double_quote = '"'

! **********************************************************************

!  process_logical_declaration() constants

! ----------------------------------------------------------------------

!  process name[=value][,name[=value]]...

character( len= *), parameter :: end_of_decl = comma // blank

! ----------------------------------------------------------------------

!  constants defining coco (integer or logical) operators, constants, etc.

! ----------------------------------------------------------------------

!  minus sign

character( len= *), parameter :: minus = '-'

!  plus sign

character( len= *), parameter :: plus = '+'

! ----------------------------------------------------------------------

!  logical binary operators

character( len= *), parameter :: or_str = '.or.'

character( len= *), parameter :: and_str = '.and.'

character( len= *), parameter :: eqv_str = '.eqv.'

character( len= *), parameter :: neqv_str = '.neqv.'

! ----------------------------------------------------------------------

!  logical uniary operator

character( len= *), parameter :: not_str = '.not.'

! ----------------------------------------------------------------------

!  logical literals

character( len= *), parameter :: true_str = '.true.'

character( len= *), parameter :: false_str = '.false.'

! ----------------------------------------------------------------------

!  the archaic versions of the relational operators

character( len= *), parameter :: dot_eq = '.eq.'

character( len= *), parameter :: dot_ne = '.ne.'

character( len= *), parameter :: dot_gt = '.gt.'

character( len= *), parameter :: dot_ge = '.ge.'

character( len= *), parameter :: dot_le = '.le.'

character( len= *), parameter :: dot_lt = '.lt.'

!  the modern versions of the relational operators

character( len= *), parameter :: ch_eq = '=='

character( len= *), parameter :: ch_ne = '/='

character( len= *), parameter :: ch_gt = '>'

character( len= *), parameter :: ch_ge = '>='

character( len= *), parameter :: ch_le = '<='

character( len= *), parameter :: ch_lt = '<'

! ----------------------------------------------------------------------

!  strings used to declare symbol names and values

! ----------------------------------------------------------------------

!  equal sign

character( len= *), parameter :: equals = '='

!  open parenthesis

character( len= *), parameter :: open_paren = '('

!  close parenthesis

character( len= *), parameter :: close_paren = ')'

! ----------------------------------------------------------------------

!  directives which must appear in the setfile

! ----------------------------------------------------------------------

!  alter directive

character( len= *), parameter :: alter_str = 'alter:'

!  directory declaration

character( len= *), parameter :: directory_str = 'directory'

!  edit directive allows changing the edit mode from the setfile

character( len= *), parameter :: edit_str = 'edit:'

!  docfile declaration

character( len= *), parameter :: docfile_str = 'docfile'

!  logfile declaration

character( len= *), parameter :: logfile_str = 'logfile'

!  number directive controls placing "! file: line" strings on source lines

character( len= *), parameter :: number_str = 'number:'

!  parens directive sets warning when actual arguments don't have enclosing parenthesis

character( len= *), parameter :: parens_str = 'parens:'

!  summary directive allows changing the summary mode from the setfile

character( len= *), parameter :: summary_str = 'summary:'

!  verbose directive allows changing the verbose mode from the setfile

character( len= *), parameter :: verbose_str = 'verbose:'

!  warn directive allows changing the warning mode from the setfile

character( len= *), parameter :: warn_str = 'warn:'

!  wrap directive allows changing the wrap length from the setfile

character( len= *), parameter :: wrap_str = 'wrap:'

! ----------------------------------------------------------------------

!  directives which may appear in the setfile or source file

! ----------------------------------------------------------------------

!  integer declaration

character( len= *), parameter :: integer_str = 'integer::'

!  integer constant declaration

character( len= *), parameter :: integer_constant_str = 'integer,parameter::'

!  logical declaration

character( len= *), parameter :: logical_str = 'logical::'

!  logical constant declaration

character( len= *), parameter :: logical_constant_str = 'logical,parameter::'

! ----------------------------------------------------------------------

!  directives which must appear in the source file

! ----------------------------------------------------------------------

!  include directive

character( len= *), parameter :: include_str = 'include'

! ----------------------------------------------------------------------

!  stop directive

character( len= *), parameter :: stop_str = 'stop'

! ----------------------------------------------------------------------

!  message directive

character( len= *), parameter :: message_str = 'message'

! ----------------------------------------------------------------------

!  if directive

character( len= *), parameter :: if_str = 'if('

! ----------------------------------------------------------------------

!  elseif directive

character( len= *), parameter :: elseif_str = 'elseif('

! ----------------------------------------------------------------------

!  )then must close an if( or elseif(

character( len= *), parameter :: then_str = ')then'

! ----------------------------------------------------------------------

!  else directive

character( len= *), parameter :: else_str = 'else'

! ----------------------------------------------------------------------

!  endif directive

character( len= *), parameter :: endif_str = 'endif'

! ----------------------------------------------------------------------

!  directives which are extensions

! ----------------------------------------------------------------------

!  ifdef directive

character( len= *), parameter :: ifdef_str = 'ifdef('

! ----------------------------------------------------------------------

!  ifndef directive

character( len= *), parameter :: ifndef_str = 'ifndef('

! ----------------------------------------------------------------------

!  undefine directive

character( len= *), parameter :: undefine_str = 'undefine::'

! ----------------------------------------------------------------------

!  macro declaration

character( len= *), parameter :: macro_str = 'macro::'

! ----------------------------------------------------------------------

!  assert directive

character( len= *), parameter :: assert_str = 'assert'

! ----------------------------------------------------------------------

!  assertif directive

character( len= *), parameter :: assertif_str = 'assertif('

! ----------------------------------------------------------------------

!  dump directive

character( len= *), parameter :: dump_str = 'dump'

!  options directive

character( len= *), parameter :: options_str = 'options'

!  report directive

character( len= *), parameter :: report_str = 'report'

! ----------------------------------------------------------------------

!  text directive

character( len= *), parameter :: text_str = 'text::'

!  copy directive

character( len= *), parameter :: copy_str = 'copy::'

!  copyif directive

character( len= *), parameter :: copyif_str = 'copyif('

! ----------------------------------------------------------------------

!  doc directive (the end doc directive is defined in process_doc_directive() )

character( len= *), parameter :: doc_str = 'doc'

!  endfile directive

character( len= *), parameter :: endfile_str = 'endfile'

!  output directive

character( len= *), parameter :: output_str = 'output'

! ----------------------------------------------------------------------

!  these strings are replaced when editing source lines

! ----------------------------------------------------------------------

!  provide file, line, date, time strings in programs

character( len= *), parameter :: file_str = '?file?'

character( len= *), parameter :: line_str = '?line?'

character( len= *), parameter :: date_str = '?date?'

character( len= *), parameter :: time_str = '?time?'

!  provide coco rcs id string in programs

character( len= *), parameter :: coco_str = '?coco?'

! ----------------------------------------------------------------------

!  on directive

character( len= *), parameter :: on_str = 'on'

!  off directive

character( len= *), parameter :: off_str = 'off'

! **********************************************************************

!  possible states encountered during execution

! ----------------------------------------------------------------------

!  codes for possible alter states

integer, parameter :: alter_none = 0

integer, parameter :: alter_delete = 1

integer, parameter :: alter_blank = 2

integer, parameter :: alter_shift_1 = 3

integer, parameter :: alter_shift_0 = 4

integer, parameter :: alter_shift_3 = 5

! ----------------------------------------------------------------------

!  wrap lengths

integer, parameter :: wrap_none = -1

integer, parameter :: wrap_off = huge( 0)

! ----------------------------------------------------------------------

!  codes for possible symbol types

integer, parameter :: type_none = 0

integer, parameter :: type_integer = 1

integer, parameter :: type_logical = 2

integer, parameter :: type_macro = 3

integer, parameter :: type_text = 4

! ----------------------------------------------------------------------

!  codes for possible if construct phases

integer, parameter :: outside_block = 0

integer, parameter :: if_block = 1

integer, parameter :: elseif_block = 2

integer, parameter :: else_block = 3

! **********************************************************************

!  communication with getopt()

! ----------------------------------------------------------------------

!  getopt() 'no more arguments'

integer, parameter :: end_of_args = -1

!  getopt() 'not in optltrs'

character( len= *), parameter :: unknown_option = '?'

! **********************************************************************

!  coco types

! **********************************************************************

!  coco files and search paths

! ----------------------------------------------------------------------

!  file type

type :: file_t

   integer :: logical_unit

   character( len= filename_len) :: name_str

   character( len= format_len) :: format_str

   character( len= free_format_len), pointer :: line

   integer :: io_status

   integer :: lines_transfered

   logical :: named_file

   logical :: create

end type file_t

! ----------------------------------------------------------------------

!  search location type

type :: path_t

   character( len= filename_len) :: name_str

   integer :: times_accessed

   type( path_t), pointer :: next

end type path_t

! **********************************************************************

!  these derived types are used to store coco constants or variables

! ----------------------------------------------------------------------

!  type stores a coco symbol & value

type :: symbol_t

   integer :: type_code

   character( len= symbol_name_len) :: name_str

   logical :: defined

   logical :: constant

   logical :: predefined

   integer :: integer_value

   logical :: logical_value

   character( len= symbol_name_len), pointer, dimension( :) :: dummy_args

   character( len= buffer_len) :: macro_value

   character( len= free_format_len), pointer, dimension( :) :: text_lines

   type( symbol_t), pointer :: next

end type symbol_t

! **********************************************************************

!  if_t stores the state of an if block

! ----------------------------------------------------------------------

!  if_t

type :: if_t

   logical :: now_selected

   logical :: ever_selected

   integer :: phase

   type( if_t), pointer :: nested

   type( if_t), pointer :: enclosing

end type if_t

! **********************************************************************

!  state_t stores a set of coco options

! ----------------------------------------------------------------------

!  state_t

type :: state_t

   integer :: alter_state

   logical :: edit_date

   logical :: edit_file

   logical :: edit_source

   logical :: edit_integers

   logical :: edit_macros

   logical :: number_source

   logical :: args_in_parens

   logical :: report_extensions

   logical :: print_summary

   logical :: warn_undeclared

   logical :: verbose_mode

   integer :: wrap_length

end type state_t

! **********************************************************************

!  report_t stores coco statistics

! ----------------------------------------------------------------------

!  report_t records the source and sink of lines

type :: report_t

   integer :: input_lines

   integer :: input_files

   integer :: include_files

   integer :: coco_lines

   integer :: selected_lines

   integer :: elided_lines

   integer :: text_blocks

   integer :: text_lines

   integer :: copied_lines

end type report_t

! **********************************************************************

!  coco variables

! **********************************************************************

!  option swtiches

! ----------------------------------------------------------------------

type( state_t) :: options
! ----------------------------------------------------------------------

!  report totals of event counts

type( report_t) :: total

! ----------------------------------------------------------------------

!  if construct outside any if construct

type( if_t), target :: outside_any_if_construct

!  if construct status

type( if_t), pointer :: if_construct

! ----------------------------------------------------------------------

!  coco symbols are stored in a singly linked list

type( symbol_t), pointer :: first_symbol

! ----------------------------------------------------------------------

!  mark when non constants are used to provide a value for a constant

logical :: all_constants

! ----------------------------------------------------------------------

!  signal when reading the setfile

logical :: processing_set_file = .true.

! **********************************************************************

!  coco file name variables

! ----------------------------------------------------------------------

!  input file, output file, or setfile

! ----------------------------------------------------------------------

!  the (first) input file

type( file_t), target :: input_file 

! ----------------------------------------------------------------------

!  the output file

type( file_t), target :: output_file

! ----------------------------------------------------------------------

!  the docfile

type( file_t), target :: doc_file

! ----------------------------------------------------------------------

!  the setfile

type( file_t), target :: set_file

! ----------------------------------------------------------------------

!  the logfile

type( file_t) :: log_file

! ----------------------------------------------------------------------

!  point to current input file for error messages

type( file_t), pointer :: current_file

! ----------------------------------------------------------------------

!  a list of source files for reports

type( file_t), allocatable, dimension(:) :: source_file_list

!  number of filenames

integer :: number_of_names

! ----------------------------------------------------------------------

!  the input/output line buffer

character( len= free_format_len), target :: line

!  the logfile line buffer

character( len= free_format_len), target :: log_line

! ----------------------------------------------------------------------

!  list of include directories is initially . only

type( path_t), pointer :: first_directory

! ----------------------------------------------------------------------

!  communicate with getopt()

! ----------------------------------------------------------------------

!  getopt() string returning non-option letter words

character( len= filename_len) :: optarg = null_string

! ----------------------------------------------------------------------

!  number of command line args

integer :: nargs

!  count command line words

integer :: optind = 0

! **********************************************************************

!  coco local

! ----------------------------------------------------------------------

!  loop index of filename args

   integer :: this_input

! **********************************************************************

!  coco text

! **********************************************************************

continue

! ----------------------------------------------------------------------

!  initialize coco program variables

   call initialize_coco()

! ----------------------------------------------------------------------

!  process command line to get options and filenames

   call process_command_line

! ----------------------------------------------------------------------

!  see if setfile exists and process it if it does

   call seek_set_file

! ----------------------------------------------------------------------

!  set option to default values if the command line or the setfile hasn't set them

   call set_option_defaults

! ----------------------------------------------------------------------

!  open the output file but link not current_file

   call open_file( output_file)

! **********************************************************************

!  read all input file(s)

   number_of_input_files: select case( number_of_names)

   case( 0, 1) number_of_input_files

!  process the input file

      call process_input_file( input_file)

!  end of input

   case default number_of_input_files

!  process several input files

      read_all_files: do this_input = 2, number_of_names

!  process the input using coco default units

         call process_input_file( source_file_list( this_input) )

!  repeat for each input file

      enddo read_all_files

!  end of input

   end select number_of_input_files

! **********************************************************************

!  if the output file has content, copy the setfile to it

   made_output: if( output_file% lines_transfered > 0 )then

!  mark the setfile in the output (whether it is present or not)

      line = mark_set_file

      call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  if processed a setfile

      append_set_file: if( set_file% named_file )then

!  copy setfile file to output

         call copy_set_file

!  if processed setfile

      endif append_set_file

   endif made_output

! ----------------------------------------------------------------------

!  if processed a docfile

   close_doc_file: if( doc_file% named_file )then

!  copy setfile file to output

      call close_file( doc_file)

!  if processed docfile

   endif close_doc_file

! ----------------------------------------------------------------------

!  close the output file

   call close_file( output_file)

! ----------------------------------------------------------------------

!  report to logfile

   summary: if( options% print_summary )then

      call write_report

   endif summary

! ----------------------------------------------------------------------

!  close the logfile

   call close_file( log_file)

! ----------------------------------------------------------------------

!  coco exit

stop 'coco normal exit'

! **********************************************************************

!  coco library

! **********************************************************************

contains

! **********************************************************************
! **********************************************************************

!  initialize_coco() prepares coco for execution

subroutine initialize_coco()

! **********************************************************************

!  initialize_coco() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  initialize derived types values

   options = state_t( alter_none, .true., .true., .true., .true., .true., &
                      .false., .false., .false., .true., .false., .false., wrap_none)

   total = report_t( 0, 0, 0, 0, 0, 0, 0, 0, 0)

!  magic if-block outside any if block

   outside_any_if_construct = if_t( .true., .true., outside_block, null(), null() )

!  files

   input_file = file_t( input_unit, '<stdin>', string_fmt, null(), 0, 0, .false., .false.)

   output_file = file_t( output_unit, '<stdout>', string_fmt, null(), 0, 0, .false., .true.)

   doc_file = file_t( doc_unit, null_string, string_fmt, null(), 0, 0, .false., .true.)

   set_file = file_t( set_unit, null_string, string_fmt, null(), 0, 0, .false., .false.)

   log_file = file_t( error_unit, '<stderr>', integer_fmt, null(), 0, 0, .false., .true.)

!  initialize pointers

   nullify( if_construct)

   if_construct => outside_any_if_construct

   nullify( first_symbol)

   nullify( current_file)

   log_file% line => log_line

   output_file% line => line

   nullify( first_directory)

! ----------------------------------------------------------------------

!  initialize_coco() exit

return

! **********************************************************************

!  initialize_coco()

end subroutine initialize_coco

! **********************************************************************
! **********************************************************************

!  %%% open and close files, write logfile messages, parse command line

! **********************************************************************
! **********************************************************************

!  open_file() open a file and remark

subroutine open_file( this_file)

! **********************************************************************

!  open_file() interface

! ----------------------------------------------------------------------

!  the file to be opened

type( file_t), target, intent( inout) :: this_file

! **********************************************************************

!  open_file() constants

! ----------------------------------------------------------------------

!  open for reading or writing

   character( len= *), parameter :: read_action = 'READ'

   character( len= *), parameter :: write_action = 'WRITE'

!  open existing file or create a new one

   character( len= *), parameter :: read_status = 'OLD'

   character( len= *), parameter :: write_status = 'REPLACE'

! **********************************************************************

!  open_file() local

! ----------------------------------------------------------------------

!  open for reading or writing

   character( len= io_specifier_len) :: open_status

   character( len= io_specifier_len) :: open_action

! **********************************************************************

!  open_file() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  open the file if file is named

   file_has_name: if( this_file% named_file )then

!  establish open parameters for reading or writing

      direction: if( this_file% create )then

         open_status = write_status

         open_action = write_action

      else direction

         open_status = read_status

         open_action = read_action

      endif direction

!  open this file

      open( unit= this_file% logical_unit, &
            file= this_file% name_str, &
            status= open_status, &
            action= open_action, &
            iostat= this_file% io_status)

      current_file => this_file

      named_status: if( this_file% io_status > 0 )then

         call msg_quit( "can't open file: " // trim( this_file% name_str) )

      elseif( options% verbose_mode )then named_status

         call msg_continue( "opened file: " // trim( this_file% name_str) )

         nullify( current_file)

      endif named_status

   endif file_has_name

!  the logfile is never the current input file, since it receives error messages

   current_input_only: select case( this_file% logical_unit)

   case( input_unit, set_unit, read_unit: )

      current_file => this_file

      this_file% line => line

   case( output_unit)

      this_file% line => line

   case( log_unit)

      this_file% line => log_line

   end select current_input_only

! ----------------------------------------------------------------------

!  open_file() exit

return

! **********************************************************************

!  open_file()

end subroutine open_file

! **********************************************************************
! **********************************************************************

!  open_scratch() open an unformatted scratch file

subroutine open_scratch( this_file)

! **********************************************************************

!  open_scratch() interface

! ----------------------------------------------------------------------

!  the scratch file to be opened

type( file_t), target, intent( inout) :: this_file

! **********************************************************************

!  open_scratch() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  open the file

   open( unit= this_file% logical_unit, &
         status= 'SCRATCH', &
         action= 'READWRITE', &
         form= 'UNFORMATTED', &
         iostat= this_file% io_status)

   scratch_status: if( this_file% io_status > 0 )then

      current_file => this_file

      call msg_quit( "can't open scratch file: <scratch>")

   endif scratch_status

!  link to line buffer

   this_file% line => line

! ----------------------------------------------------------------------

!  open_scratch() exit

return

! **********************************************************************

!  open_scratch()

end subroutine open_scratch

! **********************************************************************
! **********************************************************************

!  close_file() close a file and remark

subroutine close_file( this_file)

! **********************************************************************

!  close_file() interface

! ----------------------------------------------------------------------

!  the file to be closed

type( file_t), target, intent( inout) :: this_file

! **********************************************************************

!  close_file() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  close the named file

   close_named: if( this_file% named_file )then

      close( unit= this_file% logical_unit, &
             status= 'KEEP', &
             iostat= this_file% io_status)

      logfile_close: if( this_file% logical_unit == log_unit )then

         this_file% logical_unit = error_unit

      endif logfile_close

      current_file => this_file

      close_status: if( this_file% io_status > 0 )then

         call msg_quit( "can't close file: " // trim( this_file% name_str) )

      elseif( options% verbose_mode )then close_status

        call msg_continue( "closed file: " // trim( this_file% name_str) )

      endif close_status

   endif close_named

!  file is not connected

   nullify( current_file)

! ----------------------------------------------------------------------

!  close_file() exit

return

! **********************************************************************

!  close_file()

end subroutine close_file

! **********************************************************************
! **********************************************************************

!  close_scratch() close a file and remark

subroutine close_scratch( this_file)

! **********************************************************************

!  close_scratch() interface

! ----------------------------------------------------------------------

!  the scratch file to be closed

type( file_t), target, intent( inout) :: this_file

! **********************************************************************

!  close_scratch() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  close the scratch file

   close( unit= this_file% logical_unit, &
          status= 'DELETE', &
          iostat= this_file% io_status)

   close_status: if( this_file% io_status > 0 )then

      current_file => this_file

      call msg_quit( "can't close scratch file: <scratch>")

   endif close_status

! ----------------------------------------------------------------------

!  close_scratch() exit

return

! **********************************************************************

!  close_scratch()

end subroutine close_scratch

! **********************************************************************
! **********************************************************************

!  set_option_defaults() set options to their default values

subroutine set_option_defaults

! **********************************************************************

!  Some options are initially set to absurd values in order to allow
!  the command line option to override the corresponding setfile directive.
!  These options need to be set to useful values after the setfile
!  has been initially read but before coco further executes.

!  These options are: the alter mode and the wrap length.

!  The options selected are also made mutually consistent.

! **********************************************************************

!  set_option_defaults() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  if the command line or the setfile hasn't set the alter state, set it to the default

   alter_default: if( options% alter_state == alter_none )then

      options% alter_state = alter_shift_3

   endif alter_default

!  if the command line or the setfile hasn't set the wrap length, set it to the default

   wrap_default: if( options% wrap_length == wrap_none )then

      options% wrap_length = free_format_len

   endif wrap_default

! ----------------------------------------------------------------------

!  ensure the correct relationship among the options

! ----------------------------------------------------------------------

!  edit specific items only if editing generally

   options% edit_date = options% edit_date .and. options% edit_source
   options% edit_file = options% edit_file .and. options% edit_source
   options% edit_integers = options% edit_integers .and. options% edit_source
   options% edit_macros = options% edit_macros .and. options% edit_source

! ----------------------------------------------------------------------

!  report specific items only if reporting generally

   options% args_in_parens = options% args_in_parens .and. options% print_summary
   options% report_extensions = options% report_extensions .and. options% print_summary
   options% warn_undeclared = options% warn_undeclared .and. options% print_summary
   options% verbose_mode = options% verbose_mode .and. options% print_summary

! ----------------------------------------------------------------------

!  set_option_defaults() exit

return

! **********************************************************************

!  set_option_defaults()

end subroutine set_option_defaults

! **********************************************************************
! **********************************************************************

!  msg_quit() process error and stop

subroutine msg_quit( msg)

! **********************************************************************

!  msg_quit() interface

! ----------------------------------------------------------------------

!  the error message

character( len= *), intent( in) :: msg

! **********************************************************************

!  msg_quit() local

! ----------------------------------------------------------------------

!  strings conatining the line number and iostat of the failed operation

   character( len= conversion_len) :: number_str

   character( len= conversion_len) :: iostat_str

! **********************************************************************

!  msg_quit() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  if file is associated with this error

   file_msg: if( associated( current_file) )then

!  if a line is associated with this error

      line_msg: if( associated( current_file% line) )then

         write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( current_file% line)

      endif line_msg

!  if io error caused this error

      io_error: if( current_file% io_status > 0 )then

!  decode line number & iostat

         write( unit= number_str, fmt= conversion_fmt) current_file% lines_transfered

         write( unit= iostat_str, fmt= conversion_fmt) current_file% io_status

!  write error message with file data

         write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( 'coco error: file: ' &
                // trim( current_file% name_str) // ', line: ' // trim( adjustl( number_str)) &
                // ', ' // ', iostat: ' // trim( adjustl( iostat_str)) // ': ' // msg)

!  if io error caused not this error

      else io_error

!  decode line number

         write( unit= number_str, fmt= conversion_fmt) current_file% lines_transfered

!  write error message with file data

         write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( 'coco error: file: ' &
                // trim( current_file% name_str) // ', line: ' // trim( adjustl( number_str)) &
                // ', ' // msg)

      endif io_error

!  if file associated not with this error

   else file_msg

!  write error message without file data

      write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( 'coco error: ' // msg)

   endif file_msg

! ----------------------------------------------------------------------

!  msg_quit() exit

stop 'coco error exit'

! **********************************************************************

!  msg_quit()

end subroutine msg_quit

! **********************************************************************
! **********************************************************************

!  msg_continue() print message or continue processing

subroutine msg_continue( msg)

! **********************************************************************

!  msg_continue() interface

! ----------------------------------------------------------------------

!  the warning or informational message

character( len= *), intent( in) :: msg

! **********************************************************************

!  msg_continue() local

! ----------------------------------------------------------------------

!  string containing the current input line number

   character( len= conversion_len) :: number_str

! **********************************************************************

!  msg_continue() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  if file associated with this message

   file_msg: if( associated( current_file) )then

!  decode line number

      write( unit= number_str, fmt= conversion_fmt) current_file% lines_transfered

!  write message with file data

      write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( 'coco message: file: ' &
             // trim( current_file% name_str) // ', line: ' // trim( adjustl( number_str)) &
             // ': ' // msg)

!  if file associated not with this message

   else file_msg

!  write message without file data

      write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( 'coco message: ' // msg)

   endif file_msg

! ----------------------------------------------------------------------

!  msg_continue() exit

return

! **********************************************************************

!  msg_continue()

end subroutine msg_continue                                       

! **********************************************************************
! **********************************************************************

!  process_command_line() process command line

subroutine process_command_line

! **********************************************************************

!  process_command_line calls getopt() to get any options, then
!  process_command_line gets filenames from the command line

! **********************************************************************

!  default coco filenames constants

! **********************************************************************

!  input filename constants

! ----------------------------------------------------------------------

!  suffix used to construct input filename if one name is on the command line

   character( len= *), parameter :: input_suffix = dot // 'fpp'

! **********************************************************************

!  output filename constants

! ----------------------------------------------------------------------

!  suffix used to construct output filename if one name is on the command line

   character( len= *), parameter :: output_suffix = dot // 'f90'

! **********************************************************************

!  setfile constants

! ----------------------------------------------------------------------

!  suffix used to construct setfile name if name is on the command line

   character( len= *), parameter :: set_suffix = dot // 'set'

! **********************************************************************

!  suffix length

   integer, parameter :: suffix_len = &
         max( len( input_suffix), len( output_suffix), len( set_suffix) )

! **********************************************************************

!  other command line constants

! ----------------------------------------------------------------------

!  coco usage (error message)

   character( len= *), parameter :: usage_msg = &
      'usage: coco [ -V  | -h | [[opts] [ basename | output input [...]]]'

! **********************************************************************

!  coco communicate with getopt()

! ----------------------------------------------------------------------

!  valid option letters

   character( len= *), parameter :: opt_letters = 'a:dD:efF:hiI:l:mnprsuvVw: '

! **********************************************************************

!  process_command_line local

! ----------------------------------------------------------------------

!  getopt() option letter

   integer :: optltr

!  input filenames

   character( len= buffer_len) :: argword

!  dot divides basename and suffix

   integer :: basename_len

!  loop through input filenames

   integer :: this_word

!  allocation status

   integer :: astat

! **********************************************************************

!  process_command_line() text

! ----------------------------------------------------------------------

continue

!  get number of command line args

   nargs = command_argument_count()

!  do until end of args is returned

   optltr = getopt( opt_letters)

!  process options

   cl_options: do while( optltr /= end_of_args)

!  select which option

      which_option: select case( char( optltr))

! ----------------------------------------------------------------------

!  set the alter state

      case( 'a') which_option

         call process_alter_option( optarg)

!  turn off ?date? & ?time? editing

      case( 'd') which_option

         options% edit_date = .false.

!  declare a symbol

      case( 'D') which_option

         call process_symbol_option( optarg)

!  turn off source editing

      case( 'e') which_option

         options% edit_source = .false.

!  turn off ?file? & ?line? editing

      case( 'f') which_option

         options% edit_file = .false.

!  write a documentation file

      case( 'F') which_option

         doc_file% name_str = optarg

         doc_file% named_file = .true.

         call open_file( doc_file)

!  help

      case( 'h') which_option

         call print_help

         stop 'coco normal exit'

!  turn off ?integer? & ?logical? editing

      case( 'i') which_option

         options% edit_integers = .false.

!  set directories to search for include files

      case( 'I') which_option

         call process_include_option( optarg)

!  set logfile (NOTE: optarg has len= filename_len, so no overflow can occur.)

      case( 'l') which_option

         log_file% logical_unit = log_unit

         log_file% name_str = optarg

         log_file% named_file = .true.

         call open_file( log_file)

!  turn off ?macro? editing

      case( 'm') which_option

         options% edit_macros = .false.

!  turn on line numbers

      case( 'n') which_option

         options% number_source = .true.

!  turn on (arg) checking

      case( 'p') which_option

         options% args_in_parens = .true.

!  turn on reporting extensions

      case( 'r') which_option

         options% report_extensions = .true.

!  turn off summary report

      case( 's') which_option

         options% print_summary = .false.

!  turn off undefined report

      case( 'u') which_option

         options% warn_undeclared = .true.

!  turn on verbose

      case( 'v') which_option

         options% verbose_mode = .true.

!  print coco version data

      case( 'V') which_option

         write( unit= error_unit, fmt= string_fmt) coco_rcs_id

         stop 'coco normal exit'

!  turn off wrapping source output

      case( 'w') which_option

         call process_wrap_value( optarg)

!  command line error

      case default which_option

         write( unit= error_unit, fmt= string_fmt) usage_msg

         stop 'coco normal exit'

      end select which_option

! ----------------------------------------------------------------------

      optltr = getopt( opt_letters)

   enddo cl_options

! ----------------------------------------------------------------------

!  the rest of the command line words (if any) must be file names

! ----------------------------------------------------------------------

!  number of command line args left unprocessed

   args_left: if( optarg == unknown_option )then

      number_of_names = nargs - optind

      optind = optind + 1

      no_more_args: if( number_of_names > 0 )then

         call get_command_argument( number= optind, value= optarg)

      endif no_more_args

   else args_left

      number_of_names = nargs - optind + 1

   endif args_left

! ----------------------------------------------------------------------

!  process filenames

   filenames: select case( number_of_names)

! ----------------------------------------------------------------------

!  one filename arg

   case( 1) filenames

!  check that basename is not too long

      base_too_long: if( ( len_trim( optarg) + suffix_len) > filename_len )then

         call msg_quit( 'filename too long: ' // trim( optarg) )

      endif base_too_long

!  use basename to make input filename

      input_file% logical_unit = read_unit

      input_file% named_file = .true.

      input_file% name_str = trim( optarg) // input_suffix

!  use basename to make output filename

      output_file% logical_unit = write_unit

      output_file% named_file = .true.

      output_file% name_str = trim( optarg) // output_suffix

!  use basename to make setfile filename

      set_file% logical_unit = set_unit

      set_file% named_file = .true.

      set_file% name_str = trim( optarg) // set_suffix

! ----------------------------------------------------------------------

!  more than one filename arg

   case( 2: ) filenames

!  read source from read_unit

      input_file% logical_unit = read_unit

      input_file% named_file = .true.

!  allocate source file list

      allocate( source_file_list( 2: number_of_names), stat= astat)

      alloc_error: if( astat > 0 )then

         call msg_quit( "can't allocate input file array")

      endif alloc_error

!  check that output name is not too long

      output_too_long: if( len_trim( optarg) > filename_len )then

         call msg_quit( 'output filename too long: ' // trim( optarg) )

      endif output_too_long

!  set up output file

      output_file% logical_unit = write_unit

      output_file% named_file = .true.

      output_file% name_str = optarg

!  compute setfile name

      basename_len = index( output_file% name_str, dot, back= .true.)

      no_dot: if( basename_len == 0 )then

         basename_len = len_trim( output_file% name_str) + len( dot)

      endif no_dot

!  check that setfile name is not too long

      set_too_long: if( basename_len + suffix_len > filename_len )then

         call msg_quit( 'setfile name too long: ' // trim( output_file% name_str) // set_suffix )

      endif set_too_long

!  set up setfile

      set_file% logical_unit = set_unit

      set_file% named_file = .true.

      set_file% name_str = output_file% name_str( : basename_len - len( dot)) // set_suffix

!  record input files in source file list

      list_inputs: do this_word = 2, number_of_names

!  establish the components of this input file except the name

         source_file_list( this_word) = input_file

!  get next arg string

         optind = optind + 1

         call get_command_argument( number= optind, value= argword)

!  check that output name is not too long

         next_too_long: if( len_trim( argword) > filename_len )then

            call msg_quit( 'input filename too long: ' // trim( optarg) )

         endif next_too_long

         source_file_list( this_word)% name_str = argword

      enddo list_inputs

!  only possible values

   end select filenames

! ----------------------------------------------------------------------

!  process_command_line() exit

return

! **********************************************************************

!  process_command_line()

end subroutine process_command_line

! **********************************************************************
! **********************************************************************

!  getopt() return next known option from command line or unknown

integer function getopt( optstring)

! **********************************************************************

!  getopt() interface

! ----------------------------------------------------------------------

!  the string of valid option letters

character( len= *), intent( in) :: optstring

! **********************************************************************

!  getopt() constants

! ----------------------------------------------------------------------

!  special characters

   character( len= *), parameter :: dash = '-'

   character( len= *), parameter :: colon = ':'

! **********************************************************************

!  getopt() local

! ----------------------------------------------------------------------

!  argument buffer

   character( len= filename_len) :: optword

!  index in optstring

   integer :: index_optstring

! **********************************************************************

!  getopt() text

continue

! ----------------------------------------------------------------------

!  initialize for next option

   check_inc: if( optind >= nargs )then

      optarg = unknown_option
      getopt = end_of_args

      return

   endif check_inc

! ----------------------------------------------------------------------

!  get next option

   optind = optind + 1

   call get_command_argument( number= optind, value= optword)

!  if word is not -?

   not_an_option: if( optword( 1: 1) /= dash )then

      optarg = optword
      getopt = end_of_args

      return

!  if word is --

   elseif( optword( 2: 2) == dash )then not_an_option

      optarg = unknown_option
      getopt = end_of_args

      return

   endif not_an_option

! ----------------------------------------------------------------------

!  optword is -x (not --)

   index_optstring = index( optstring, optword( 2: 2))

   is_opt: if( index_optstring > 0 )then

!  if this optltr must have another word

      opt_string: if( optstring( index_optstring + 1: index_optstring + 1) == colon )then

!  it can be separated by a blank

         next_word: if( optword( 3: 3) == blank )then

            optind = optind + 1
            call get_command_argument( number= optind, value= optarg)

!  or not be separated by a blank

         else next_word

            optarg = optword( 3: )

         endif next_word

      endif opt_string

      getopt = ichar( optword( 2: 2))

!  if this optltr must not have another word

   else is_opt

      optarg = optword
      getopt = ichar( unknown_option)

   endif is_opt

! ----------------------------------------------------------------------

!  getopt() exit

return

! **********************************************************************

!  getopt()

end function getopt                                               

! **********************************************************************
! **********************************************************************

!  %%% process particular command line options

! **********************************************************************
! **********************************************************************

!  process_alter_option() process alter arguments

subroutine process_alter_option( alter_opt)

! **********************************************************************

!  process_alter_option() interface

! ----------------------------------------------------------------------

!  the alter option from the command line

character( len= *), intent( in) :: alter_opt                        

! **********************************************************************

!  entry: alter_opt is command line arg following -a
!         "d" | "b" | "0" | "1" | "3"

!  exit: alter_opt is processed or error exit

! **********************************************************************

!  process_alter_option() constants

! ----------------------------------------------------------------------

!  possible alter option strings

   character( len= *), parameter :: delete_str = 'd'

   character( len= *), parameter :: blank_str = 'b'

   character( len= *), parameter :: shift0_str = '0'

   character( len= *), parameter :: shift1_str = '1'

   character( len= *), parameter :: shift3_str = '3'

! **********************************************************************

!  process_alter_option() local

! ----------------------------------------------------------------------

!  decoding the option is done in lower case which may require a case change

   character( len= 1) :: lower_case_opt

! **********************************************************************

!  process_alter_option() text

continue

! ----------------------------------------------------------------------

!  check for unknown option

   too_long: if( len_trim( alter_opt) > 1 )then

      call msg_quit( "unknown -a option: " // trim( alter_opt) )

   endif too_long

!  force arg to lower case

   fix_case: select case( alter_opt( 1: 1))

   case( 'A': 'Z') fix_case

      lower_case_opt = achar( iachar( alter_opt( 1: 1)) + change_case)

   case default fix_case

      lower_case_opt = alter_opt

   end select fix_case

! ----------------------------------------------------------------------

!  legal alter argument or error

   alter_value_str: select case( lower_case_opt)

!  alter delete

   case( delete_str) alter_value_str

      options% alter_state = alter_delete

!  alter blank

   case( blank_str) alter_value_str

      options% alter_state = alter_blank

!  alter shift1

   case( shift1_str) alter_value_str

      options% alter_state = alter_shift_1

!  alter shift0

   case( shift0_str) alter_value_str

      options% alter_state = alter_shift_0

!  alter shift3

   case( shift3_str) alter_value_str

      options% alter_state = alter_shift_3

!  unknown alter code ( not one of { b, d, 0, 1, 3 } )

   case default alter_value_str

      call msg_quit( "unknown -a option: " // trim( alter_opt) )

!  legal alter statement or error

   end select alter_value_str

! ----------------------------------------------------------------------

!  process_alter_option() exit

return

! **********************************************************************

!  process_alter_option()

end subroutine process_alter_option

! **********************************************************************
! **********************************************************************

!  process_symbol_option() process define arguments

subroutine process_symbol_option( symbol_opt)

! **********************************************************************

!  process_symbol_option() interface

! ----------------------------------------------------------------------

!  the symbol string from the command line

character( len= *), intent( in) :: symbol_opt                       

! **********************************************************************

!  entry: symbol_opt is string following -D { log | int=val }

!  exit: symbol_opt is processed or error exit

! **********************************************************************

!  process_symbol_option() constants

! ----------------------------------------------------------------------

   character( len= *), parameter :: log_value_str = '=.true.'

! **********************************************************************

!  process_symbol_option() local

! ----------------------------------------------------------------------

!  find characters

   integer :: next_char

!  construct a declaration string to process

   character( len= filename_len) :: decl_str

! **********************************************************************

!  process_symbol_option() text

continue

! ----------------------------------------------------------------------

!  force names to lower case

   each_char: do next_char = 1, len( symbol_opt)

      to_lower: select case( symbol_opt( next_char: next_char))

      case( 'A': 'Z') to_lower

         decl_str( next_char: next_char) = achar( iachar( symbol_opt( next_char: next_char)) + change_case)

      case default to_lower

         decl_str( next_char: next_char) = symbol_opt( next_char: next_char)

      end select to_lower

   enddo each_char

! ----------------------------------------------------------------------

!  an equal sign must separate a value from the name

   next_char = index( decl_str, equals)

!  if there's an equals, it's an integer

   int_or_log: if( next_char > 0 )then

!  declare the integer constant

      call process_integer_constant( decl_str)

!  if there's no equals, it's a logical ( = .true.)

   else int_or_log

!  construct the logical (it's true)

      decl_str = trim( decl_str) // log_value_str

!  declare the logical constant

      call process_logical_constant( decl_str)

!  integer or logical

   endif int_or_log

! ----------------------------------------------------------------------

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "defined symbol from command line: " // trim( symbol_opt))

   endif extensions

! ----------------------------------------------------------------------

!  process_symbol_option() exit

return

! **********************************************************************

!  process_symbol_option()

end subroutine process_symbol_option

! **********************************************************************
! **********************************************************************

!  print_help() write summary report to specified unit

subroutine print_help

! **********************************************************************

!  entry: in response to -h command line option

!  exit: print help message

! **********************************************************************

!  print_help() constants

! ----------------------------------------------------------------------

!  the help message

   character( len= *), dimension( 20), parameter :: help_msg = (/ &
                      ' -a ?           set alter state, ? = { b, d, 0, 1, 3}       ', &
                      ' -d             turn off ?date? & ?time? editing            ', &
                      ' -D name[=n]    declare integer or logical constant         ', &
                      ' -e             turn off all source editing                 ', &
                      ' -f             turn off ?file? & ?line? editing            ', &
                      ' -F file        write documentation text to file            ', &
                      ' -h             print this help message and quit            ', &
                      ' -i             turn off ?integer? & ?logical? editing      ', &
                      ' -I dir         search dir for include files (after .)      ', &
                      ' -l file        write log messages to file (default stderr) ', &
                      ' -m             turn off ?macro? editing                    ', &
                      ' -n             print line numbers on source lines          ', &
                      ' -p             warn when actual args might need parenthesis', &
                      ' -r             report all use of extensions                ', &
                      ' -s             work silently, print no summary             ', &
                      ' -u             report setfile symbols undefined in source  ', &
                      ' -v             report file opening and closing             ', &
                      ' -V             print coco version and quit                 ', &
                      ' -w n           wrap lines to n columns (0= off, 72= fixed) ', &
                      ' --             optionally separate options from file names ' /)

! **********************************************************************

!  print_help() local

! ----------------------------------------------------------------------

!  implied do variable

   integer :: do_idx

! **********************************************************************

!  print_help() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

   write( unit= error_unit, fmt= string_fmt) ( trim( help_msg( do_idx)), do_idx= 1, size( help_msg))

! ----------------------------------------------------------------------

!  print_help() exit

return

! **********************************************************************

!  print_help()

end subroutine print_help

! **********************************************************************
! **********************************************************************

!  process_include_option() process include directory options

subroutine process_include_option( directory_opt)

! **********************************************************************

!  process_include_option() interface

! ----------------------------------------------------------------------

!  the directory string from the command line

character( len= *), intent( in) :: directory_opt                   

! **********************************************************************

!  entry: directory_opt is a directory to be added to the list
!         of directories to be searched for inlcude files

!  exit: directory_opt is on the list

! **********************************************************************

!  process_include_option() local

! ----------------------------------------------------------------------

!  point to a directory type

   type( path_t), pointer :: path_ptr

!  unquote the directory name if needed

   character( len= filename_len) :: directory_str

!  lengths of quoted string

   integer :: quoted_len

!  lengths of unquoted string

   integer :: unquoted_len

! **********************************************************************

!  process_include_option() text

continue

! ----------------------------------------------------------------------

!  if the directory is quoted, unquote it

   is_quoted: select case( directory_opt( 1: 1) )

   case( single_quote, double_quote) is_quoted

      call unquote_string( directory_opt, directory_str, quoted_len, unquoted_len)

      null_unquoted: if( quoted_len == 0 .or. unquoted_len == 0 )then

         call msg_quit( "null name passed to -I option")

      endif null_unquoted

   case default is_quoted

      directory_str = directory_opt

   end select is_quoted

! ----------------------------------------------------------------------

!  if name is already on the path

   nullify( path_ptr)

   call seek_directory( directory_str, path_ptr)

   on_list_or_add: if( associated( path_ptr) )then

      call msg_continue( "redundant include directory ignored: " // trim( directory_opt) )

   else on_list_or_add

      call add_directory( directory_str)

!  if reporting use of extensions

      extensions: if( options% report_extensions )then

         call msg_continue( "added include directory from command line: " // trim( directory_opt) )

      endif extensions

   endif on_list_or_add

! ----------------------------------------------------------------------

!  process_include_option() exit

return

! **********************************************************************

!  process_include_option()

end subroutine process_include_option

! **********************************************************************
! **********************************************************************

!  write_options() write summary report to specified unit

subroutine write_options

! **********************************************************************

!  write_options() constants

! ----------------------------------------------------------------------

!  possible alter states

integer, parameter :: lower_alter = min( alter_delete, alter_blank, alter_shift_1, alter_shift_0, alter_shift_3)

integer, parameter :: upper_alter = max( alter_delete, alter_blank, alter_shift_1, alter_shift_0, alter_shift_3)

!  possible alter state labels

character( len= 16), dimension( lower_alter: upper_alter), parameter :: alter_labels = (/ &
                                       'deleted         ', &
                                       'blank line      ', &
                                       'initial !       ', &
                                       'shifted 1 + !   ', &
                                       'shifted 3 + !?> ' /)

! **********************************************************************

!  write_options() local

! ----------------------------------------------------------------------

!  construct output lines

   character( len= source_line_len) :: output_line

! **********************************************************************

!  write_options() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  write a header

   write( unit= log_file% logical_unit, fmt= log_file% format_str) "coco options:"

! ----------------------------------------------------------------------

!  identify the alter state

   check_index: select case( options% alter_state)

   case( lower_alter: upper_alter) check_index

      output_line = 'alter state causes lines to be ' // alter_labels( options% alter_state)

   case default check_index

      output_line = 'alter state is undefined'

   end select check_index

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether editing date & time

   edit_date_time: if( options% edit_date )then

      output_line = 'editing ' // date_str // ' and ' // time_str // ' strings'

   else edit_date_time

      output_line = 'not editing ' // date_str // ' and ' // time_str // ' strings'

   endif edit_date_time

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether editing file & line

   edit_file_line: if( options% edit_file )then

      output_line = 'editing ' // file_str // ' and ' // line_str // ' strings'

   else edit_file_line

      output_line = 'not editing ' // file_str // ' and ' // line_str // ' strings'

   endif edit_file_line

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether editing at all

   edit_control: if( options% edit_source )then

      output_line = 'editing source lines'

   else edit_control

      output_line = 'not editing source lines'

   endif edit_control

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether editing integers and logicals

   edit_ints_logs: if( options% edit_integers )then

      output_line = 'editing integer and logicals'

   else edit_ints_logs

      output_line = 'not editing integer and logicals'

   endif edit_ints_logs

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether editing macros

   edit_macro: if( options% edit_macros )then

      output_line = 'editing macros'

   else edit_macro

      output_line = 'not editing macros'

   endif edit_macro

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether warning whether args have enclosing parens

   edit_parens: if( options% args_in_parens )then

      output_line = 'warning when macro & text actual args are not enclosed in parenthesis'

   else edit_parens

      output_line = 'not warning when macro & text actual args are not enclosed in parenthesis'

   endif edit_parens

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether reporting extensions

   rpt_exten: if( options% report_extensions )then

      output_line = 'reporting coco extensions'

   else rpt_exten

      output_line = 'not reporting coco extensions'

   endif rpt_exten

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether printing coco summary

   rpt_prt: if( options% print_summary )then

      output_line = 'printing coco report'

   else rpt_prt

      output_line = 'not printing coco report'

   endif rpt_prt

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether printing "! file: line" on source lines

   rpt_number: if( options% number_source )then

      output_line = 'numbering source lines'

   else rpt_number

      output_line = 'not numbering source lines'

   endif rpt_number

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether warning when symbol is in setfile but not source

   rpt_warn: if( options% warn_undeclared )then

      output_line = 'warning when setfile symbols not in source'

   else rpt_warn

      output_line = 'not warning when setfile symbols not in source'

   endif rpt_warn

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether verbose mode is on

   rpt_verbose: if( options% verbose_mode )then

      output_line = 'verbose mode is on'

   else rpt_verbose

      output_line = 'verbose mode is off'

   endif rpt_verbose

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line)

!  identify whether verbose mode is on

   rpt_wrap: if( options% wrap_length /= wrap_off )then

      output_line = 'wrapping source lines at length '

   else rpt_wrap

      output_line = 'not wrapping source lines'

   endif rpt_wrap

   write( unit= log_file% logical_unit, fmt= log_file% format_str) trim( output_line) // blank, options% wrap_length

! ----------------------------------------------------------------------

!  write_options() exit

return

! **********************************************************************

!  write_options()

end subroutine write_options

! **********************************************************************
! **********************************************************************

!  write_report() write summary report to specified unit

subroutine write_report

! **********************************************************************

!  write_report() local

! ----------------------------------------------------------------------

!  print date and time in header

   character( len= 8) :: today

   character( len= 10) :: now

! ----------------------------------------------------------------------

!  print include path

   type( path_t), pointer :: path_ptr

!  search lists for symbols not defined in the coco program proper

   type( symbol_t), pointer :: symbol_ptr

! ----------------------------------------------------------------------

!  print input files

   integer :: this_file

! **********************************************************************

!  write_report() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  banner includes the date and time

   call date_and_time( date= today, time= now)

   write( unit= log_file% logical_unit, fmt= log_file% format_str) 'coco executed: ' // today // blank // now

! ----------------------------------------------------------------------

!  identify the setfile

   write( unit= log_file% logical_unit, fmt= log_file% format_str) 'setfile: ' // trim( set_file% name_str)

! ----------------------------------------------------------------------

!  identify the output file

   write( unit= log_file% logical_unit, fmt= log_file% format_str) 'output: ' // trim( output_file% name_str)

! ----------------------------------------------------------------------

!  identify the input file(s)

   one_or_more: if( allocated( source_file_list) )then

      write( unit= log_file% logical_unit, fmt= log_file% format_str, advance= 'NO') 'input:'

      more_than_one: do this_file = 2, number_of_names

         write( unit= log_file% logical_unit, fmt= log_file% format_str, advance= 'NO') &
                blank // trim( source_file_list( this_file)% name_str)

      enddo more_than_one

      write( unit= log_file% logical_unit, fmt= log_file% format_str)

   else one_or_more

      write( unit= log_file% logical_unit, fmt= log_file% format_str) 'input: ' // trim( input_file% name_str)

   endif one_or_more

! ----------------------------------------------------------------------

!  identify the document file, if there is one

   got_doc: if( doc_file% named_file )then

      write( unit= log_file% logical_unit, fmt= log_file% format_str) &
         'docfile: ' // trim( doc_file% name_str), doc_file% lines_transfered

   endif got_doc

! ----------------------------------------------------------------------

!  identify the include path

   write( unit= log_file% logical_unit, fmt= log_file% format_str, advance= 'NO') 'include path: .'

   nullify( path_ptr)

   path_ptr => first_directory

   inc_path: do while( associated( path_ptr) )

      write( unit= log_file% logical_unit, fmt= directory_fmt, advance= 'NO') &
             blank // trim( path_ptr% name_str) // open_paren, path_ptr% times_accessed, close_paren

      path_ptr => path_ptr% next

   enddo inc_path

!  end line using null string

   write( unit= log_file% logical_unit, fmt= log_file% format_str)

! ----------------------------------------------------------------------

!  if undefined symbols report is requested

   undefined_complaints: if( options% warn_undeclared )then

!  complain about any integers or logicals declared in the setfile but not in source

      nullify( symbol_ptr)

      symbol_ptr => first_symbol

      search_syms: do while( associated( symbol_ptr))

         found_sym: if( symbol_ptr% predefined )then

            write( unit= log_file% logical_unit, fmt= log_file% format_str) &
               'symbol declared in setfile but not in any source file: ' // trim( symbol_ptr% name_str)

         endif found_sym

         symbol_ptr => symbol_ptr% next

      enddo search_syms

   endif undefined_complaints

! ----------------------------------------------------------------------

!  number of files read

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'files read: ', total% input_files

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'include files read: ', total% include_files

!  number of setfile lines read

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'set lines read: ', set_file% lines_transfered

!  number of coco lines read

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'coco lines read: ', total% coco_lines

!  number of source lines read

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'source lines read: ', total% input_lines

!  number of lines written

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'source lines written: ', output_file% lines_transfered

!  number of selected lines written

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'selected source lines: ', total% selected_lines

!  number of elided lines

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'elided source lines: ', total% elided_lines

!  number of text blocks read

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'text blocks read: ', total% text_blocks

!  number of text lines read

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'text lines read: ', total% text_lines

!  number of text lines written

   write( unit= log_file% logical_unit, fmt= log_file% format_str) &
          'text lines written: ', total% copied_lines

! ----------------------------------------------------------------------

!  write_report() exit

return

! **********************************************************************

!  write_report()

end subroutine write_report

! **********************************************************************
! **********************************************************************

!  process_wrap_value() set the wrap length option length

subroutine process_wrap_value( number_str)

! **********************************************************************

!  process_wrap_value() interface

! ----------------------------------------------------------------------

!  the wrap length string from the command line

character( len= *), intent( in) :: number_str

! **********************************************************************

!  process_wrap_value() local

! ----------------------------------------------------------------------

!  index of a non-digit

   integer :: char_idx

!  convert string to characters

   character( len= conversion_len) :: conversion_str

!  test proposed wrap length

   integer :: wrap_len

! **********************************************************************

!  process_wrap_value() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  legal wrap string value or error

   bad_char: if( char_idx > 0 )then

      char_idx = verify( number_str, digit_chars)

      call msg_quit( "bad character in wrap length: " // trim( number_str))

   endif bad_char

! ----------------------------------------------------------------------

!  convert chacacters to integer

   conversion_str = number_str

   conversion_str = adjustr( conversion_str)

   read( unit= conversion_str, fmt= conversion_fmt) wrap_len

! ----------------------------------------------------------------------

!  check proposed wrap length

   set_within_bounds: select case( wrap_len)

!  do no wrapping

   case( 0) set_within_bounds

      options% wrap_length = wrap_off

!  within bounds so accept as if

   case( fixed_format_len: free_format_len) set_within_bounds

      options% wrap_length = wrap_len

!  out of bounds so set within bounds

   case default set_within_bounds

      call msg_continue( "invalid wrap length set within bounds: " // trim( number_str))

      options% wrap_length = min( free_format_len, max( fixed_format_len, wrap_len))

   end select set_within_bounds

! ----------------------------------------------------------------------

!  process_wrap_value() exit

return

! **********************************************************************

!  process_wrap_value()

end subroutine process_wrap_value

! **********************************************************************
! **********************************************************************

!  %%% seek and process the setfile (if any)

! **********************************************************************
! **********************************************************************

!  seek_setfile() write summary report to specified unit

subroutine seek_set_file                                            

! **********************************************************************

!  seek_set_file constants

! ----------------------------------------------------------------------

!  default set_file name

   character( len= *), parameter :: default_name = 'coco.set'

! **********************************************************************

!  seek_set_file() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  see if there is name for the set_file

   named_set_file: if( set_file% named_file )then

!  inquire by file looking for named setfile

      inquire( file= set_file% name_str, exist= set_file% named_file, iostat= set_file% io_status)

      inq_named: if( set_file% io_status > 0 )then

         call msg_quit( "can't inquire setfile: " // trim( set_file% name_str))

      endif inq_named

!  done checking setfile name

   endif named_set_file

! ----------------------------------------------------------------------

!  no named setfile so try to find the default setfile

   default_set_file: if( .not. set_file% named_file )then

!  inquire by file looking for default setfile

      inquire( file= default_name, exist= set_file% named_file, iostat= set_file% io_status)

      inq_default: if( set_file% io_status > 0 )then

         call msg_quit( "can't inquire default setfile: " // default_name)

      endif inq_default

!  if found the default setfile ensure the variable correctly specifies it

      use_default: if( set_file% named_file )then

         set_file% logical_unit = set_unit

         set_file% name_str = default_name

      endif use_default

   endif default_set_file

! ----------------------------------------------------------------------

!  if have setfile, open it, process it, close it

   read_set_file: if( set_file% named_file )then

      call process_set_file

   else read_set_file

      set_file% name_str = '<no setfile>'

   endif read_set_file

   processing_set_file = .false.

! ----------------------------------------------------------------------

!  seek_set_file() exit

return

! **********************************************************************

!  seek_set_file()

end subroutine seek_set_file

! **********************************************************************
! **********************************************************************

!  process_sefile() open, process, close the coco setfile

subroutine process_set_file                                         

! **********************************************************************

!  process_set_file() steps

!  1. open the setfile
!  2. open the set scratch file
!  2. read the setfile line by line
!  3. call blank_compress_lower_case() to construct a coco statement
!  4. ignore coco comments
!  5. call process_set_statement() to process coco set statement
!  6. close setfile

! **********************************************************************

!  process_set_file() local

! ----------------------------------------------------------------------

!  process the setfile statement by statement

   character( len= buffer_len) :: set_statement

! ----------------------------------------------------------------------

!  signal complete statement

   logical :: complete

! **********************************************************************

!  process_set_file() text

continue

! ----------------------------------------------------------------------

!  open the setfile for reading

   call open_file( set_file)

! ----------------------------------------------------------------------

!  count files processed

   total% input_files = total% input_files + 1

!  as if finished a complete statement at beginning of file

   complete = .true.

! ----------------------------------------------------------------------

!  main read setfile lines loop

   read_lines: do

! ----------------------------------------------------------------------

!  read a setfile line

      read( unit= set_file% logical_unit, fmt= set_file% format_str, iostat= set_file% io_status) set_file% line

      read_set: if( set_file% io_status > 0 )then

         call msg_quit( "can't read setfile: " // trim( set_file% name_str))

      endif read_set

! ----------------------------------------------------------------------

!  read until end of file

      read_eof: if( set_file% io_status < 0 )then

!  reset statement processing for the next file

         call blank_compress_lower_case( set_statement, null_string)

!  if in a statement continuation sequence

         premature_eof: if( .not. complete )then

            call msg_quit( "end of file encountered within a continuation sequence")

         endif premature_eof

!  exit the read lines loop

         exit read_lines

      endif read_eof

!  count setfile lines

      set_file% lines_transfered = set_file% lines_transfered + 1

! ----------------------------------------------------------------------

!  process setfile lines or error if source lines

      coco_line: if( line( : len( coco_key)) == coco_key )then

!  count coco lines

         total% coco_lines = total% coco_lines + 1

!  process setfile lines, ignore coco comments

         coco_statement: if( is_coco_statement( line( len( coco_key) + 1: )) )then

! ----------------------------------------------------------------------

!  read a complete statement line by line

            call gather_coco_statement( line, set_statement, complete)

!  if not yet a complete statement go get the rest of it

            get_statement: if( .not. complete )then

               cycle read_lines

            endif get_statement

!  process the complete setfile statement

            call process_set_statement( set_statement)

!  process setfile lines, ignore coco comments

         endif coco_statement

!  source line in setfile

      else coco_line

         call msg_quit( "source lines are not allowed in the setfile")

!  end processing set statements

      endif coco_line

! ----------------------------------------------------------------------

!  end main read setfile lines loop

   enddo read_lines

   total% input_lines = total% input_lines + set_file% lines_transfered

! ----------------------------------------------------------------------

!  close the setfile

   call close_file( set_file)

! ----------------------------------------------------------------------

!  process_set_file() exit

return

! **********************************************************************

!  process_set_file()

end subroutine process_set_file

! **********************************************************************
! **********************************************************************

!  copy_set_file() copy setfile to output file

subroutine copy_set_file                                             

! **********************************************************************

!  copy_set_file() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  open the set file

   call open_file( set_file)

! ----------------------------------------------------------------------

!  copy each line

   copy_lines: do

!  read a line

      read( unit= set_file% logical_unit, fmt= set_file% format_str, &
            iostat= set_file% io_status) set_file% line

      read_set_file: if( set_file% io_status > 0 )then

         call msg_quit( "can't copy setfile")

      endif read_set_file

!  read entire scratch file

      set_eof: if( set_file% io_status < 0 )then

         exit copy_lines

      endif set_eof

!  write a line

      call write_coco_line( output_file)

   enddo copy_lines

! ----------------------------------------------------------------------

!  close the setfile

   call close_file( set_file)

! ----------------------------------------------------------------------

!  copy_set_file() exit

return

! **********************************************************************

!  copy_set_file()

end subroutine copy_set_file

! **********************************************************************
! **********************************************************************

!  %%% process statements many of which may appear only in the setfile

! **********************************************************************
! **********************************************************************

!  process_set_statement() process set line

subroutine process_set_statement( set_stmt)

! **********************************************************************

!  process_set_statement() interface

! ----------------------------------------------------------------------

!  the statement string from the set file

character( len= *), intent( in) :: set_stmt                       

! **********************************************************************

!  entry: set_stmt is blank_compress_lower_case set statement past the coco key
!         "alter:..." | "integer..." | "logical..." | "directory'...'" |
!         "wrap:..." | "edit:..." | "warn:..." | "logfile'...'" |
!         "summary:..." | "verbose:..." | "doc'...'" | "parens:..."

!  exit: set_stmt is processed or error exit

! **********************************************************************

!  process_set_statement() text

continue

! ----------------------------------------------------------------------

!  catergorize setfile statement: alter, integer, logical, directory, wrap

! ----------------------------------------------------------------------

!  if the directive is an alter directive

   which_directive: if( set_stmt( : len( alter_str)) == alter_str )then

      call process_alter_directive( set_stmt( len( alter_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile integer declaration

   elseif( set_stmt( : len( integer_str)) == integer_str )then which_directive

      call process_integer_declaration( set_stmt( len( integer_str) + 1: ) )

! ----------------------------------------------------------------------

!  integer constant declaration

   elseif( set_stmt( : len( integer_constant_str)) == integer_constant_str )then which_directive

      call process_integer_constant( set_stmt( len( integer_constant_str) + 1: ))

! ----------------------------------------------------------------------

!  if the directive is a setfile logical declaration

   elseif( set_stmt( : len( logical_str)) == logical_str )then which_directive

      call process_logical_declaration( set_stmt( len( logical_str) + 1: ) )

! ----------------------------------------------------------------------

!  logical constant declaration

   elseif( set_stmt( : len( logical_constant_str)) == logical_constant_str )then which_directive

      call process_logical_constant( set_stmt( len( logical_constant_str) + 1: ))

! ----------------------------------------------------------------------

!  if the directive is a setfile directory directive

   elseif( set_stmt( : len( directory_str)) == directory_str )then which_directive

      call process_directory_directive( set_stmt( len( directory_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile docfile directive

   elseif( set_stmt( : len( docfile_str)) == docfile_str )then which_directive

      call process_docfile_directive( set_stmt( len( docfile_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile logfile directive

   elseif( set_stmt( : len( logfile_str)) == logfile_str )then which_directive

      call process_logfile_directive( set_stmt( len( logfile_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile wrap directive

   elseif( set_stmt( : len( wrap_str)) == wrap_str )then which_directive

      call process_wrap_directive( set_stmt( len( wrap_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile edit directive

   elseif( set_stmt( : len( edit_str)) == edit_str )then which_directive

      call process_edit_directive( set_stmt( len( edit_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile summary directive

   elseif( set_stmt( : len( summary_str)) == summary_str )then which_directive

      call process_summary_directive( set_stmt( len( summary_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile number: directive

   elseif( set_stmt( : len( number_str)) == number_str )then which_directive

      call process_number_directive( set_stmt( len( number_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile parens directive

   elseif( set_stmt( : len( parens_str)) == parens_str )then which_directive

      call process_paren_directive( set_stmt( len( parens_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile warn directive

   elseif( set_stmt( : len( warn_str)) == warn_str )then which_directive

      call process_warn_directive( set_stmt( len( warn_str) + 1: ) )

! ----------------------------------------------------------------------

!  if the directive is a setfile verbose directive

   elseif( set_stmt( : len( verbose_str)) == verbose_str )then which_directive

      call process_verbose_directive( set_stmt( len( verbose_str) + 1: ) )

! ----------------------------------------------------------------------

!  otherwise complain about the unknown directive

   else which_directive

      call msg_quit( "unknown setfile directive: " // trim( set_stmt) )

!  catergorize setfile statement: alter or integer or logical

   endif which_directive

! ----------------------------------------------------------------------

!  process_set_statement() exit

return

! **********************************************************************

!  process_set_statement()

end subroutine process_set_statement

! **********************************************************************
! **********************************************************************

!  process_alter_directive() process alter directives

subroutine process_alter_directive( alter_dir)

! **********************************************************************

!  process_alter_directive() interface

! ----------------------------------------------------------------------

!  the alter directive from the setfile

character( len= *), intent( in) :: alter_dir                        

! **********************************************************************

!  entry: alter_dir is blank_compress_lower_case alter directive past the colon
!         "delete" | "blank" | "shift0" | "shift1" | "shift3"

!  exit: alter_dir is processed or error exit

! **********************************************************************

!  process_alter_directive() constants

! ----------------------------------------------------------------------

!  possible alter strings

   character( len= *), parameter :: delete_str = 'delete'

   character( len= *), parameter :: blank_str = 'blank'

   character( len= *), parameter :: shift0_str = 'shift0'

   character( len= *), parameter :: shift1_str = 'shift1'

   character( len= *), parameter :: shift3_str = 'shift3'

! **********************************************************************

!  process_alter_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_alter_statements = .false.

! **********************************************************************

!  process_alter_directive() text

continue

! ----------------------------------------------------------------------

!  only one alter directive per setfile

   too_many_alters: if( too_many_alter_statements )then

      call msg_quit( "too many alter statements")

   else too_many_alters

      too_many_alter_statements = .true.

   endif too_many_alters

!  if the alter state has not been set from the command line

   not_set: if( options% alter_state == alter_none )then

! ----------------------------------------------------------------------

!  legal alter statement or error

! ----------------------------------------------------------------------

!  decode alter state

      alter_value_str: select case( alter_dir)

 !  alter delete

     case( delete_str) alter_value_str

         options% alter_state = alter_delete

!  alter blank

      case( blank_str) alter_value_str

         options% alter_state = alter_blank

!  alter shift1

      case( shift1_str) alter_value_str

         options% alter_state = alter_shift_1

!  alter shift0

      case( shift0_str) alter_value_str

         options% alter_state = alter_shift_0

!  alter shift3

      case( shift3_str) alter_value_str

         options% alter_state = alter_shift_3

!  unknown alter

      case default alter_value_str

         call msg_quit( "unknown alter directive: " // trim( alter_dir))

!  legal alter statement or error

      end select alter_value_str

   endif not_set

! ----------------------------------------------------------------------

!  process_alter_directive() exit

return

! **********************************************************************

!  process_alter_directive()

end subroutine process_alter_directive

! **********************************************************************
! **********************************************************************

!  process_directory_directive() process include directory options

subroutine process_directory_directive( directory_dir)

! **********************************************************************

!  process_directory_directive() interface

! ----------------------------------------------------------------------

!  the directory directive from the setfile

character( len= *), intent( in) :: directory_dir                    

! **********************************************************************

!  entry: directory_opt is a directory to be added to the list
!         of directories to be searched for inlcude files

!  exit: directory_opt is on the list

! **********************************************************************

!  process_directory_directive() local

! ----------------------------------------------------------------------

!  point to a directory type

   type( path_t), pointer :: directory_ptr

!  the name of a directory

   character( len= filename_len) :: name_str

!  count length of quoted string

   integer :: directive_len

!  count length of unquoted string

   integer :: name_len

! **********************************************************************

!  process_directory_directive() text

continue

! ----------------------------------------------------------------------

!  unquote string to find path string

   call unquote_string( directory_dir, name_str, directive_len, name_len )

   no_name_str: if( name_len == 0 .or. directive_len == 0 )then

      call msg_quit( "no directory name: " // trim( directory_dir) )

   endif no_name_str

!  verify no extra characters beyond name

   extra_chars: if( directory_dir( directive_len + 1: ) /= blank )then

      call msg_quit( "extra characters after directory path name: " // trim( directory_dir))

   endif extra_chars

! ----------------------------------------------------------------------

!  if name is already on the path

   call seek_directory( name_str, directory_ptr)

   on_list_or_add: if( associated( directory_ptr) )then

      call msg_continue( "redundant include directory ignored: " // trim( directory_dir) )

!  if name is not already on the path

   else on_list_or_add

      call add_directory( name_str)

!  if reporting use of extensions

      extensions: if( options% report_extensions )then

         call msg_continue( "added include directory from setfile: " // trim( directory_dir) )

      endif extensions

   endif on_list_or_add

! ----------------------------------------------------------------------

!  process_directory_directive() exit

return

! **********************************************************************

!  process_directory_directive()

end subroutine process_directory_directive

! **********************************************************************
! **********************************************************************

!  seek_directory() return a pointer to directory_str or null()

subroutine seek_directory( name_str, directory_ptr)

! **********************************************************************

!  seek_directory() interface

! ----------------------------------------------------------------------

!  the name of the directory to seek

character( len= *), intent( in) :: name_str

!  a pointer to the directory entry if found or null()

type( path_t), pointer :: directory_ptr

! **********************************************************************

!  entry: directory_str is a directory to be added to the list
!         of directories to be searched for inlcude files

!  exit: directory_str is on the list

! **********************************************************************

!  seek_directory() text

continue

! ----------------------------------------------------------------------

!  search from beginning to end of path list

   nullify( directory_ptr)

   directory_ptr => first_directory

!  if the name is already in the path

   scan_path: do while( associated( directory_ptr) )

      found_name: if( name_str == directory_ptr% name_str )then

         exit scan_path

      endif found_name

      directory_ptr => directory_ptr% next

   enddo scan_path

! ----------------------------------------------------------------------

!  seek_directory() exit

return

! **********************************************************************

!  seek_directory()

end subroutine seek_directory

! **********************************************************************
! **********************************************************************

!  add_directory() return a pointer to directory_str or null()

subroutine add_directory( directory_str)

! **********************************************************************

!  add_directory() interface

! ----------------------------------------------------------------------

!  the name of the directory to add to the directory list

character( len= *), intent( in) :: directory_str                    

! **********************************************************************

!  entry: directory_str is a directory to be added to the list
!         of directories to be searched for inlcude files

!  exit: directory_str is on the list

! **********************************************************************

!  add_directory() local

! ----------------------------------------------------------------------

!  end of linked list, null() if no linked list yet

   type( path_t), save, pointer :: current_directory => null()

!  check allocation status

   integer :: astat

! **********************************************************************

!  add_directory() text

continue

! ----------------------------------------------------------------------

!  append to list

   start_or_append: if( associated( first_directory) )then

      allocate( current_directory% next, stat= astat)

      append_status: if( astat > 0 )then

         call msg_quit( "can't append to path list: " // trim( directory_str) )

      endif append_status

      current_directory => current_directory% next

!  start list

   else start_or_append

      allocate( first_directory, stat= astat)

      start_status: if( astat > 0 )then

         call msg_quit( "can't start path list: " // trim( directory_str) )

      endif start_status

      current_directory => first_directory

   endif start_or_append

!  update new entry

   current_directory% name_str = directory_str

   current_directory% times_accessed = 0

   nullify( current_directory% next)

! ----------------------------------------------------------------------

!  add_directory() exit

return

! **********************************************************************

!  add_directory()

end subroutine add_directory

! **********************************************************************
! **********************************************************************

!  process_docfile_directive() process include docfile options

subroutine process_docfile_directive( docfile_dir)

! **********************************************************************

!  process_docfile_directive() interface

! ----------------------------------------------------------------------

!  the docfile directive from the setfile

character( len= *), intent( in) :: docfile_dir                      

! **********************************************************************

!  entry: docfile_opt is a docfile to be added to the list
!         of directories to be searched for inlcude files

!  exit: docfile_opt is on the list

! **********************************************************************

!  process_docfile_directive() local

! ----------------------------------------------------------------------

!  the name of the file to be opened

   character( len= filename_len) :: docfile_name

!  the length of the quoted string

   integer :: quoted_len

!  the length of the unquoted string

   integer :: unquoted_len

!  count number of some statements to disallow more than one

   logical, save :: too_many_docfile_statements = .false.

! **********************************************************************

!  process_docfile_directive() text

continue

! ----------------------------------------------------------------------

!  only one docfile statement per setfile

   too_many_docfiles: if( too_many_docfile_statements )then

      call msg_quit( "too many docfile statements")

   else too_many_docfiles

      too_many_docfile_statements = .true.

   endif too_many_docfiles

!  unquote string on directive

   call unquote_string( docfile_dir, docfile_name, unquoted_len, quoted_len)

   no_name: if( unquoted_len == 0 .or. quoted_len == 0 )then

      call msg_quit( "no name found on docfile directive: " // trim( docfile_dir) )

   endif no_name

!  verify no extra characters beyond name

   extra_chars: if( docfile_dir( unquoted_len + 1: ) /= blank )then

      call msg_quit( "extra characters after docfile file name: " // trim( docfile_dir))

   endif extra_chars

! ----------------------------------------------------------------------

!  if docfile named on command line ignore the directive

   already_named: if( doc_file% named_file )then

      call msg_continue( "command line overrides setfile, docfile directive ignored: " // trim( docfile_dir) )

!  if docfile not named on command line open the named file

   else already_named

      doc_file% name_str = docfile_name

      doc_file% named_file = .true.

      call open_file( doc_file)

   endif already_named

! ----------------------------------------------------------------------

!  process_docfile_directive() exit

return

! **********************************************************************

!  process_docfile_directive()

end subroutine process_docfile_directive

! **********************************************************************
! **********************************************************************

!  process_logfile_directive() process include logfile options

subroutine process_logfile_directive( logfile_dir)

! **********************************************************************

!  process_logfile_directive() interface

! ----------------------------------------------------------------------

!  the logfile directive from the setfile

character( len= *), intent( in) :: logfile_dir                      

! **********************************************************************

!  entry: logfile_opt is a logfile to be added to the list
!         of directories to be searched for inlcude files

!  exit: logfile_opt is on the list

! **********************************************************************

!  process_logfile_directive() local

! ----------------------------------------------------------------------

!  the name of the file to be opened

   character( len= filename_len) :: logfile_name

!  the length of the quoted string

   integer :: quoted_len

!  the length of the unquoted string

   integer :: unquoted_len

!  count number of some statements to disallow more than one

   logical, save :: too_many_logfile_statements = .false.

! **********************************************************************

!  process_logfile_directive() text

continue

! ----------------------------------------------------------------------

!  only one logfile statement per setfile

   too_many_logfiles: if( too_many_logfile_statements )then

      call msg_quit( "too many logfile statements")

   else too_many_logfiles

      too_many_logfile_statements = .true.

   endif too_many_logfiles

!  unquote string on directive

   call unquote_string( logfile_dir, logfile_name, unquoted_len, quoted_len)

   no_name: if( unquoted_len == 0 .or. quoted_len == 0 )then

      call msg_quit( "no name found on logfile directive: " // trim( logfile_dir) )

   endif no_name

!  verify no extra characters beyond name

   extra_chars: if( logfile_dir( unquoted_len + 1: ) /= blank )then

      call msg_quit( "extra characters after logfile file name: " // trim( logfile_dir))

   endif extra_chars

! ----------------------------------------------------------------------

!  if logfile named on command line ignore the directive

   already_named: if( log_file% named_file )then

      call msg_continue( "command line overrides setfile, logfile directive ignored: " // trim( logfile_dir) )

!  if logfile not named on command line open the named file

   else already_named

      log_file% logical_unit = log_unit

      log_file% name_str = logfile_name

      log_file% named_file = .true.

      call open_file( log_file)

   endif already_named

! ----------------------------------------------------------------------

!  process_logfile_directive() exit

return

! **********************************************************************

!  process_logfile_directive()

end subroutine process_logfile_directive

! **********************************************************************
! **********************************************************************

!  process_wrap_directive() process wrap directives

subroutine process_wrap_directive( wrap_dir)

! **********************************************************************

!  process_wrap_directive() interface

! ----------------------------------------------------------------------

!  the wrap directive from the setfile

character( len= *), intent( in) :: wrap_dir                         

! **********************************************************************

!  entry: wrap_dir is blank_compress_lower_case wrap directive
!         it must be a number string

!  exit: wrap_dir is processed or error exit

! **********************************************************************

!  process_wrap_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_wrap_statements = .false.

! **********************************************************************

!  process_wrap_directive() text

continue

! ----------------------------------------------------------------------

!  only one wrap statement per setfile

   too_many_wraps: if( too_many_wrap_statements )then

      call msg_quit( "too many wrap statements")

   else too_many_wraps

      too_many_wrap_statements = .true.

   endif too_many_wraps

! ----------------------------------------------------------------------

!  process wrap value if not already set on command line

   set_on_cl: if( options% wrap_length == wrap_none )then

      call process_wrap_value( wrap_dir)

   endif set_on_cl

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed wrap directive from setfile: " // trim( wrap_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_wrap_directive() exit

return

! **********************************************************************

!  process_wrap_directive()

end subroutine process_wrap_directive

! **********************************************************************
! **********************************************************************

!  process_edit_directive() process edit directives

subroutine process_edit_directive( edit_dir)

! **********************************************************************

!  process_edit_directive() interface

! ----------------------------------------------------------------------

!  the edit directive from the setfile

character( len= *), intent( in) :: edit_dir

! **********************************************************************

!  entry: edit_dir is blank_compress_lower_case edit directive
!         it must be "on" or "off"

!  exit: edit_dir is processed or error exit

! **********************************************************************

!  process_edit_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_edit_statements = .false.

! **********************************************************************

!  process_edit_directive() text

continue

! ----------------------------------------------------------------------

!  only one edit statement per setfile

   too_many_edits: if( too_many_edit_statements )then

      call msg_quit( "too many edit statements")

   else too_many_edits

      too_many_edit_statements = .true.

   endif too_many_edits

! ----------------------------------------------------------------------

!  process edit switch

   on_off: if( edit_dir == on_str )then

      options% edit_source = .true.

   elseif( edit_dir == off_str)then on_off

      options% edit_source = .false.

   else on_off

      call msg_quit( "unknown option on edit directive: " // trim( edit_dir) )

   endif on_off

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed edit directive from setfile: " // trim( edit_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_edit_directive() exit

return

! **********************************************************************

!  process_edit_directive()

end subroutine process_edit_directive

! **********************************************************************
! **********************************************************************

!  process_number_directive() process number directives

subroutine process_number_directive( number_dir)

! **********************************************************************

!  process_number_directive() interface

! ----------------------------------------------------------------------

!  the number directive from the setfile

character( len= *), intent( in) :: number_dir

! **********************************************************************

!  entry: number_dir is blank_compress_lower_case number directive
!         it must be "on" or "off"

!  exit: number_dir is processed or error exit

! **********************************************************************

!  process_number_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_number_statements = .false.

! **********************************************************************

!  process_number_directive() text

continue

! ----------------------------------------------------------------------

!  only one number statement per setfile

   too_many_number: if( too_many_number_statements )then

      call msg_quit( "too many number statements")

   else too_many_number

      too_many_number_statements = .true.

   endif too_many_number

! ----------------------------------------------------------------------

!  process number switch

   on_off: if( number_dir == on_str )then

      options% number_source = .true.

   elseif( number_dir == off_str )then on_off

      options% number_source = .false.

   else on_off

      call msg_quit( "unknown option on number directive: " // trim( number_dir) )

   endif on_off

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed number directive from setfile: " // trim( number_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_number_directive() exit

return

! **********************************************************************

!  process_number_directive()

end subroutine process_number_directive

! **********************************************************************
! **********************************************************************

!  process_paren_directive() process paren directives

subroutine process_paren_directive( paren_dir)

! **********************************************************************

!  process_paren_directive() interface

! ----------------------------------------------------------------------

!  the paren directive from the setfile

character( len= *), intent( in) :: paren_dir

! **********************************************************************

!  entry: paren_dir is blank_compress_lower_case paren directive
!         it must be "on" or "off"

!  exit: paren_dir is processed or error exit

! **********************************************************************

!  process_paren_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_paren_statements = .false.

! **********************************************************************

!  process_paren_directive() text

continue

! ----------------------------------------------------------------------

!  only one paren statement per setfile

   too_many_paren: if( too_many_paren_statements )then

      call msg_quit( "too many paren statements")

   else too_many_paren

      too_many_paren_statements = .true.

   endif too_many_paren

! ----------------------------------------------------------------------

!  process paren switch

   on_off: if( paren_dir == on_str )then

      options% args_in_parens = .true.

   elseif( paren_dir == off_str )then on_off

      options% args_in_parens = .false.

   else on_off

      call msg_quit( "unknown option on paren directive: " // trim( paren_dir) )

   endif on_off

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed paren directive from setfile: " // trim( paren_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_paren_directive() exit

return

! **********************************************************************

!  process_paren_directive()

end subroutine process_paren_directive

! **********************************************************************
! **********************************************************************

!  process_summary_directive() process summary directives

subroutine process_summary_directive( summary_dir)

! **********************************************************************

!  process_summary_directive() interface

! ----------------------------------------------------------------------

!  the srap directive from the setfile

character( len= *), intent( in) :: summary_dir

! **********************************************************************

!  entry: summary_dir is blank_compress_lower_case summary directive
!         it must be "on" or "off"

!  exit: summary_dir is processed or error exit

! **********************************************************************

!  process_summary_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_summary_statements = .false.

! **********************************************************************

!  process_summary_directive() text

continue

! ----------------------------------------------------------------------

!  only one summary statement per setfile

   too_many_summarys: if( too_many_summary_statements )then

      call msg_quit( "too many summary statements")

   else too_many_summarys

      too_many_summary_statements = .true.

   endif too_many_summarys

! ----------------------------------------------------------------------

!  process summary switch

   on_off: if( summary_dir == on_str )then

      options% print_summary = .true.

   elseif( summary_dir == off_str)then on_off

      options% print_summary = .false.

   else on_off

      call msg_quit( "unknown option on summary directive: " // trim( summary_dir) )

   endif on_off

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed summary directive from setfile: " // trim( summary_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_summary_directive() exit

return

! **********************************************************************

!  process_summary_directive()

end subroutine process_summary_directive

! **********************************************************************
! **********************************************************************

!  process_verbose_directive() process verbose directives

subroutine process_verbose_directive( verbose_dir)

! **********************************************************************

!  process_verbose_directive() interface

! ----------------------------------------------------------------------

!  the verbose directive from the setfile

character( len= *), intent( in) :: verbose_dir

! **********************************************************************

!  entry: verbose_dir is blank_compress_lower_case verbose directive
!         it must be "on" or "off"

!  exit: verbose_dir is processed or error exit

! **********************************************************************

!  process_verbose_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_verbose_statements = .false.

! **********************************************************************

!  process_verbose_directive() text

continue

! ----------------------------------------------------------------------

!  only one verbose statement per setfile

   too_many_verboses: if( too_many_verbose_statements )then

      call msg_quit( "too many verbose statements")

   else too_many_verboses

      too_many_verbose_statements = .true.

   endif too_many_verboses

! ----------------------------------------------------------------------

!  process verbose switch if not already set on command line

   on_off: if( verbose_dir == on_str )then

      options% verbose_mode = .true.

   elseif( verbose_dir == off_str)then on_off

      options% verbose_mode = .false.

   else on_off

      call msg_quit( "unknown option on verbose directive: " // trim( verbose_dir) )

   endif on_off

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed verbose directive from setfile: " // trim( verbose_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_verbose_directive() exit

return

! **********************************************************************

!  process_verbose_directive()

end subroutine process_verbose_directive

! **********************************************************************
! **********************************************************************

!  process_warn_directive() process warn directives

subroutine process_warn_directive( warn_dir)

! **********************************************************************

!  process_warn_directive() interface

! ----------------------------------------------------------------------

!  the warn directive from the setfile

character( len= *), intent( in) :: warn_dir

! **********************************************************************

!  entry: warn_dir is blank_compress_lower_case warn directive
!         it must be a "on" or "off"

!  exit: warn_dir is processed or error exit

! **********************************************************************

!  process_warn_directive() local

! ----------------------------------------------------------------------

!  count number of some statements to disallow more than one

   logical, save :: too_many_warn_statements = .false.

! **********************************************************************

!  process_warn_directive() text

continue

! ----------------------------------------------------------------------

!  only one warn statement per setfile

   too_many_warns: if( too_many_warn_statements )then

      call msg_quit( "too many warn statements")

   else too_many_warns

      too_many_warn_statements = .true.

   endif too_many_warns

! ----------------------------------------------------------------------

!  process warn switch

   on_off: if( warn_dir == on_str )then

      options% warn_undeclared = .true.

   elseif( warn_dir == off_str)then on_off

      options% warn_undeclared = .false.

   else on_off

      call msg_quit( "unknown option on warn directive: " // trim( warn_dir) )

   endif on_off

!  if reporting use of extensions

   extensions: if( options% report_extensions )then

      call msg_continue( "processed warn directive from setfile: " // trim( warn_dir) )

   endif extensions

! ----------------------------------------------------------------------

!  process_warn_directive() exit

return

! **********************************************************************

!  process_warn_directive()

end subroutine process_warn_directive

! **********************************************************************
! **********************************************************************

!  %%% read and write files, and process coco lines and statements

! **********************************************************************
! **********************************************************************

!  process_input_file() reads a coco source file, recurse upon include files

recursive subroutine process_input_file( this_file)

! **********************************************************************

!  process_input_file() interface

! ----------------------------------------------------------------------

!  the file to be processed

type( file_t), target, intent( inout) :: this_file

! **********************************************************************

!  entry: source file to be processed

!  exit: source file has been processed or error

! **********************************************************************

!  process_input_file() local

! ----------------------------------------------------------------------

!  statement buffer

   character( len= buffer_len) :: statement

!  signal complete statement

   logical :: complete

! **********************************************************************

!  process_input_file() text

continue

!  open source file

   call open_file( this_file)

!  count files

   total% input_files = total% input_files + 1

!  as if finished a complete statement at beginning of file

   complete = .true.

! ----------------------------------------------------------------------

!  main read lines loop

   read_lines: do

!  read from input file

      read( unit= this_file% logical_unit, fmt= this_file% format_str, iostat= this_file% io_status) &
            this_file% line

      read_source: if( this_file% io_status > 0 )then

         call msg_quit( "can't read input file: " // trim( this_file% name_str))

      endif read_source

! ----------------------------------------------------------------------

!  read until end of file or complete statement

      read_eof: if( this_file% io_status < 0 )then

         call blank_compress_lower_case( statement, null_string)

         continuation_eof: if( .not. complete )then

            call msg_quit( "end of file encountered in a continued coco statement")

         endif continuation_eof

         exit read_lines

      endif read_eof

! ----------------------------------------------------------------------

!  count lines

      this_file% lines_transfered = this_file% lines_transfered + 1

!  process coco lines or source lines

      process_coco: if( line( : len( coco_key)) == coco_key )then

!  count lines

         total% coco_lines = total% coco_lines + 1

!  write coco line to the output

         call write_coco_line( output_file)

!  if line is not a coco comment

         coco_statement: if( is_coco_statement( line( len( coco_key) + 1: )) )then

!  read a complete statement

            call gather_coco_statement( line, statement, complete)

!  if not yet a complete statement go get the rest of it

            get_statement: if( .not. complete )then

               cycle read_lines

            endif get_statement

! ----------------------------------------------------------------------

!  process coco directives

! ----------------------------------------------------------------------

!  directive is a coco include directive

            process_directive: if( statement( : len( include_str)) == include_str )then

!  process (possibly recursive) include directives (include 'a' --> include 'b' &c)

               call process_include_directive( statement( len( include_str) + 1: ))

!  directive is a coco endfile directive

            elseif( statement == endfile_str )then process_directive

!  if line is active read no more from this input file

               active_endfile: if( if_construct% now_selected )then

                  exit read_lines

               endif active_endfile

!  process any other directive

            else process_directive

!  process other (not possibly recursive) coco directives

               call process_coco_statement( statement)

            endif process_directive

!  reset current file

            current_file => this_file

!  if line is not a coco comment

         endif coco_statement

! ----------------------------------------------------------------------

!  process source lines

      else process_coco

!  error if a source line is mixed into a continued coco statement

         continuation_error: if( .not. complete )then

            call msg_quit( "source line encountered in a continued coco statement")

         endif continuation_error

!  if within the active block of a coco if construct

         active_source: if( if_construct% now_selected )then

!  editing source lines is enabled

            edit_source: if( options% edit_source )then

!  if ? present, edit source line

               edit_line: if( index( line, arg_key) > 0 )then

                  call edit_source_line( line)

               endif edit_line

            endif edit_source

         endif active_source

!  copy source lines

         call write_source_line( output_file)

!  end processing coco lines

      endif process_coco

!  end main read lines loop

   enddo read_lines

   total% input_lines = total% input_lines + this_file% lines_transfered

! ----------------------------------------------------------------------

!  end of file

   call close_file( this_file)

! ----------------------------------------------------------------------

!  process_input_file() exit

return

! **********************************************************************

!  process_input_file()

end subroutine process_input_file

! **********************************************************************
! **********************************************************************

!  gather_coco_statement() examine lines and signal a complete statement

subroutine gather_coco_statement( line, statement, complete)

! **********************************************************************

!  gather_coco_statement() interface

! ----------------------------------------------------------------------

!  the current input file

character( len= *), intent( in) :: line

!  the statement as it is built

character( len= *), intent( inout) :: statement

!  true when a complete statement has been seen

logical, intent( out) :: complete                                   

! **********************************************************************

!  entry: statement is a line
!         "..."

!  exit: statement is accumulated, complete is true when whole

! **********************************************************************

!  gather_coco_statement() local

! ----------------------------------------------------------------------

!  count continuation lines

   integer, save :: continuation_lines = 0

! ----------------------------------------------------------------------

!  number of characters processed so far

   integer, save :: statement_len = 0

! **********************************************************************

!  gather_coco_statement() text

continue

! ----------------------------------------------------------------------

!  blank compress lower case

   call blank_compress_lower_case( statement, line( len( coco_key) + 1: ) )

!  if statement length hasn't changed and statement is not a comment

   null_stmt: if( statement_len == len_trim( statement) )then

      call msg_quit( "null statement encountered in continuation sequence")

   endif null_stmt

!  if not a complete statement yet

   statement_len = len_trim( statement)

!  last character is continuation means more to read to get a complete statement

   incomplete: if( statement( statement_len: statement_len) == continuation )then

!  if too many continuation lines

      too_many: if( continuation_lines > max_continuations )then

         call msg_quit( "too many continuations")

      endif too_many

!  count continuation lines

      continuation_lines = continuation_lines + 1

!  go get the rest of the statement

      complete = .false.

      return

   endif incomplete

!  coco statement is complete

   continuation_lines = 0

   statement_len = 0

   complete = .true.

! ----------------------------------------------------------------------

!  gather_coco_statement() exit

return

! **********************************************************************

!  gather_coco_statement()

end subroutine gather_coco_statement

! **********************************************************************
! **********************************************************************

!  is_coco_statement() true if coco statement is a coco construct and not a coco comment

logical function is_coco_statement( coco_stmt)

! **********************************************************************

!  is_coco_statement() interface

! ----------------------------------------------------------------------

!  the coco statement to be categorized

character( len= *), intent( in) :: coco_stmt                       

! **********************************************************************

!  entry: coco_stmt is coco statement past the coco key
!         "..."

!  exit: true if statement is a coco construct,
!        false if statement is a coco comment

! **********************************************************************

!  is_coco_statement() local

! ----------------------------------------------------------------------

!  locations of specific characters

   integer :: char_idx

! **********************************************************************

!  is_coco_statement() text

continue

! ----------------------------------------------------------------------

!  scan from first character (past coco key) for non-whitespace

   char_idx = verify( coco_stmt, white_space)

!  if found other than whitespace

   white_space_or_comment: if( char_idx > 0 )then

!  is construct if first character past the white space is not comment

      is_coco_statement = coco_stmt( char_idx: char_idx) /= comment

!  all whitespace

   else white_space_or_comment

!  is not a construct

      is_coco_statement = .false.

   endif white_space_or_comment

! ----------------------------------------------------------------------

!  is_coco_statement() exit

return

! **********************************************************************

!  is_coco_statement()

end function is_coco_statement

! **********************************************************************
! **********************************************************************

!  write_coco_line() write a coco line of output

subroutine write_coco_line( this_file)

! **********************************************************************

!  write_coco_line() interface

! ----------------------------------------------------------------------

!  the file to receive the output

type( file_t), intent( inout) :: this_file                        

! **********************************************************************

!  entry: out_unit is the logical unit connected to the output file
!         coco_line is a line to be written as per the current alter state

!  exit: line has been written or error exit

! **********************************************************************

!  write_coco_line() constants

! ----------------------------------------------------------------------

!  line prefix when alter state is shift 3

   character( len= *), parameter :: shift3_prefix = '!?>'

! **********************************************************************

!  write_coco_line() local

! ----------------------------------------------------------------------

!  maximum line length

   character( len= source_line_len) :: line_buffer

!  detect lines longer than free_format_len

   logical :: long_lines

! **********************************************************************

!  write_coco_line() text

continue

! ----------------------------------------------------------------------

!  write output as per alter state

   long_lines = .false.

! ----------------------------------------------------------------------

   alter_print: select case( options% alter_state)

!  delete the line

   case( alter_delete) alter_print

      return

!  blank line

   case( alter_blank) alter_print

      write( unit= this_file% logical_unit, fmt= this_file% format_str, iostat= this_file% io_status) &
             null_string

!  comment the line

   case( alter_shift_0) alter_print

      line_buffer = comment // this_file% line( len( comment) + 1: )

      write( unit= this_file% logical_unit, fmt= this_file% format_str, iostat= this_file% io_status) &
             trim( line_buffer)

!  shift one or comment the line

   case( alter_shift_1) alter_print

      line_buffer = comment // this_file% line

      long_lines = len_trim( line_buffer) > free_format_len

      write( unit= this_file% logical_unit, fmt= this_file% format_str, iostat= this_file% io_status) &
             trim( line_buffer)

!  shift three or comment the line

   case( alter_shift_3) alter_print

      line_buffer = shift3_prefix // this_file% line

      long_lines = len_trim( line_buffer) > free_format_len

      write( unit= this_file% logical_unit, fmt= this_file% format_str, iostat= this_file% io_status) &
             trim( line_buffer)

   end select alter_print

! ----------------------------------------------------------------------

!  check write iostat

   write_this_file: if( this_file% io_status > 0 )then

      call msg_quit( "error writing source output: " // trim( this_file% name_str) )

   endif write_this_file

! ----------------------------------------------------------------------

!  complain if source line too long

   len_warning: if( long_lines )then

      call msg_continue( "line longer than 132 characters:")

      call msg_continue( trim( line_buffer))

   endif len_warning

! ----------------------------------------------------------------------

!  write_coco_line() exit

return

! **********************************************************************

!  write_coco_line()

end subroutine write_coco_line

! **********************************************************************
! **********************************************************************

!  write_source_line() write a line of output

subroutine write_source_line( this_file)

! **********************************************************************

!  write_source_line() interface

! ----------------------------------------------------------------------

!  the file to receive the output

type( file_t), intent( inout) :: this_file

! **********************************************************************

!  entry: out_unit is the logical unit connected to the output file
!         source_line is the line of Fortran source to be written

!  exit: the line is written or error exit

! **********************************************************************

!  write_source_line() constants

! ----------------------------------------------------------------------

!  column to start line number

   integer, parameter :: number_len = 75

! **********************************************************************

!  write_source_line() local

! ----------------------------------------------------------------------

!  character line number

   character( len= conversion_len) :: number_str

! **********************************************************************

!  write_source_line() text

continue

! ----------------------------------------------------------------------

!  if currently printing output

   process_line: if( if_construct% now_selected )then

!  check whether to number the source line

      want_numbers: if( options% number_source )then

!  only number lines if there's nothing there now

         can_number: if( this_file% line( number_len: ) == blank )then

            write( unit= number_str, fmt= conversion_fmt) current_file% lines_transfered

            this_file% line( number_len: ) = '! ' // trim( current_file% name_str) // ': ' // adjustl( number_str)

         endif can_number

      endif want_numbers

!  write source output

      write( unit= this_file% logical_unit, fmt= this_file% format_str, iostat= this_file% io_status) &
             trim( this_file% line)

!  check for write error

      write_error: if( this_file% io_status > 0 )then

         call msg_quit( "error writing source output: " // trim( this_file% name_str) )

      endif write_error

!  count lines written

      this_file% lines_transfered = this_file% lines_transfered + 1

      total% selected_lines = total% selected_lines + 1

   else process_line

!  otherwise print as per the alter state

      call write_coco_line( this_file)

      total% elided_lines = total% elided_lines + 1

   endif process_line

! ----------------------------------------------------------------------

!  write_source_line() exit

return

! **********************************************************************

!  write_source_line()

end subroutine write_source_line

! **********************************************************************
! **********************************************************************

!  process_include_directive() process an include directive

recursive subroutine process_include_directive( include_dir)

! **********************************************************************

!  process_include_directive() interface

! ----------------------------------------------------------------------

!  the include directive

character( len= *), intent( in) :: include_dir                      

! **********************************************************************

!  entry: inc_name is inlcude file name

!  exit: inc_name is inlcude file name with directory prepended

! **********************************************************************

!  process_include_directive() constants

! ----------------------------------------------------------------------

!  mark the beginning and end of include files as per the standard

   character( len= *), parameter :: begin_inc = '??! INCLUDE '

   character( len= *), parameter :: end_inc = '??! END INCLUDE '

! **********************************************************************

!  process_include_directive() local

! ----------------------------------------------------------------------

!  file variable of file named on the include directive

   type( file_t), target :: include_file

! ----------------------------------------------------------------------

!  length of quoted include file name

   integer :: construct_len

!  length of unquoted include file name

   integer :: name_len

! **********************************************************************

!  process_include_directive() text

continue

! ----------------------------------------------------------------------

!  check the include syntax: unquote the include filename

   call unquote_string( include_dir, include_file% name_str, construct_len, name_len )

   no_name_str: if( name_len == 0 .or. construct_len == 0 )then

      call msg_quit( "no include file name: " // trim( include_dir))

   endif no_name_str

! ----------------------------------------------------------------------

!  if active block, process include directive

   active_inc: if( if_construct% now_selected )then

!  see if the include file exists

      inquire( file= include_file% name_str, &
               exist= include_file% named_file, iostat= include_file% io_status)

      inquire_error: if( include_file% io_status > 0 )then

         call msg_quit( "can't inquire include file: " // trim( include_file% name_str) )

      endif inquire_error

!  if not found, check directories

      seek_inc: if( .not. include_file% named_file )then

         call seek_include_file( include_file)

      endif seek_inc

!  if still not found, complain and quit

      no_name: if( .not. include_file% named_file )then

         call msg_quit( "can't find include file: " // trim( include_file% name_str) )

      endif no_name

! ----------------------------------------------------------------------

!  build include_file to pass to process_input_file()

! ----------------------------------------------------------------------

!  get new unit

      include_unit: if( current_file% logical_unit == input_unit )then

         include_file% logical_unit = read_unit

      else include_unit

         include_file% logical_unit = current_file% logical_unit + 1

      endif include_unit

!  include file components

      include_file% format_str = current_file% format_str

      include_file% line => null()

      include_file% io_status = 0

      include_file% lines_transfered = 0

      include_file% named_file = .true.

      include_file% create = .false.

! ----------------------------------------------------------------------

!  mark include file in output

      line = begin_inc // include_file% name_str

      call write_coco_line( output_file)

!  prepare to process include file

      total% include_files = total% include_files + 1

!  process include file

      call process_input_file( include_file)

!  relink pointers

      output_file% line => line

!  mark include file in output

      line = end_inc // include_file% name_str

      call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  if active block, process include directive

   endif active_inc

!  end processing include statement

! ----------------------------------------------------------------------

!  process_include_directive() exit

return

! **********************************************************************

!  process_include_directive()

end subroutine process_include_directive

! **********************************************************************
! **********************************************************************

!  seek_include_file() seek inlcude file in directories

subroutine seek_include_file( include_file)

! **********************************************************************

!  seek_include_file() interface

! ----------------------------------------------------------------------

!  the include file name to be sought

type( file_t), intent( inout) :: include_file

! **********************************************************************

!  entry: include_file is inlcude file name

!  exit: include_file is inlcude file name with directory prepended

! **********************************************************************

!  seek_include_file() local

! ----------------------------------------------------------------------

!  pointer to directories on path

   type( path_t), pointer :: directory

!  construct path/names to check for existance

   character( len= filename_len) :: trial_name

! **********************************************************************

!  seek_include_file() text

continue

! ----------------------------------------------------------------------

!  search list for directory/file

   nullify( directory)

   directory => first_directory

!  last directory in path is not associated

   search_list: do while( associated( directory) )

!  construct full name <directory-name><file-name>

      trial_name = trim( directory% name_str) // include_file% name_str

!  check file existence

      inquire( file= trial_name, &
               exist= include_file% named_file, iostat= include_file% io_status)

      inquire_error: if( include_file% io_status > 0 )then

         call msg_quit( "can't inquire include file: " // trim( trial_name))

      endif inquire_error

!  found file name

      name_match: if( include_file% named_file )then

!  found a file in this directory

         directory% times_accessed = directory% times_accessed + 1

!  rewrite include file name to include directory

         include_file% name_str = trial_name

         extensions: if( options% report_extensions )then

            call msg_continue( "found include file in directory: " // trim( trial_name))

         endif extensions

         exit search_list

      endif name_match

!  filename not yet found so try next directory

      directory => directory% next

   enddo search_list

! ----------------------------------------------------------------------

!  seek_include_file() exit

return

! **********************************************************************

!  seek_include_file()

end subroutine seek_include_file

! **********************************************************************
! **********************************************************************

!  %%% edit Fortran source lines

! **********************************************************************
! **********************************************************************

!  edit_source_line() edit source lines

subroutine edit_source_line( source_line)

! **********************************************************************

!  edit_source_line() interface

! ----------------------------------------------------------------------

!  source line to be edited

character( len= *), intent( inout) :: source_line

! **********************************************************************

!  entry: line is a line of Fortran source
!         with (possibly) ?file?, ?line?, ?date?, ?time?, ?integer?, ?logical?, ?macro?

!  exit: line has any ?macro? etc. strings replaced with their values

! **********************************************************************

!  edit_source_line() local

! ----------------------------------------------------------------------

!  copy of line since editing may expand the line beyond its length

   character( len= buffer_len) :: edit_line

!  make lower case copy of line

   character( len= buffer_len) :: lower_case_line

   integer :: this_char

! **********************************************************************

!  edit_source_line() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  the line may be expanded by editing, so edit a long buffer

   edit_line = source_line

!  make lower case copy of line

   scan_line: do this_char = 1, len( line)

      to_lower: select case( edit_line( this_char: this_char))

      case( 'A': 'Z') to_lower

         lower_case_line( this_char: this_char) = achar( iachar( edit_line( this_char: this_char)) &
                                                         + change_case)

      case default to_lower

         lower_case_line( this_char: this_char) = edit_line( this_char: this_char)

      end select to_lower

   enddo scan_line

! ----------------------------------------------------------------------

!  process ?coco?

   call edit_coco_strings( edit_line, lower_case_line)

! ----------------------------------------------------------------------

!  if processing ?file? & ?line?

   edit_file_and_line: if( options% edit_file )then

      call edit_file_line_strings( edit_line, lower_case_line)

   endif edit_file_and_line

! ----------------------------------------------------------------------

!  if processing ?date? & ?time?

   edit_date_and_time: if( options% edit_date )then

      call edit_date_time_strings( edit_line, lower_case_line)

   endif edit_date_and_time

! ----------------------------------------------------------------------

!  replace ?name? with the current integer value of name

   edit_integers: if( options% edit_integers )then

      call edit_integer_strings( edit_line, lower_case_line)

      call edit_logical_strings( edit_line, lower_case_line)

   endif edit_integers

! ----------------------------------------------------------------------

!  replace ?name? with the current string value of name

   edit_macros: if( options% edit_macros )then

      call edit_macro_strings( edit_line, lower_case_line)

   endif edit_macros

! ----------------------------------------------------------------------

!  remove any line length overflow

   wrap_lines: if( options% wrap_length /= wrap_off )then

      call wrap_source_line( edit_line)

   endif wrap_lines

   source_line = edit_line

! ----------------------------------------------------------------------

!  edit_source_line() exit

return

! **********************************************************************

!  edit_source_line()

end subroutine edit_source_line

! **********************************************************************
! **********************************************************************

!  edit_coco_strings() process ?coco? strings

subroutine edit_coco_strings( edit_line, lower_case_line)

! **********************************************************************

!  edit_coco_strings() interface

! ----------------------------------------------------------------------

!  the source line to be edited

character( len= *), intent( inout) :: edit_line

!  the source line in lower case to enable searches

character( len= *), intent( inout) :: lower_case_line

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?integer

!  exit: line has any ?name strings replaced with their values

! **********************************************************************

!  edit_coco_strings() local

! ----------------------------------------------------------------------

!  find substring

   integer :: search_idx

! **********************************************************************

!  edit_coco_strings() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  replace ?coco? with the coco rcs id

   search_idx = index( lower_case_line, coco_str)

   go_coco: if( search_idx > 0 )then

      call replace_substring( edit_line, lower_case_line, coco_str, coco_rcs_id, search_idx)

      remark_coco: if( options% report_extensions )then

         call msg_continue( "edited coco: " // coco_rcs_id)

      endif remark_coco

   endif go_coco

! ----------------------------------------------------------------------

!  edit_coco_strings() exit

return

! **********************************************************************

!  edit_coco_strings()

end subroutine edit_coco_strings

! **********************************************************************
! **********************************************************************

!  edit_file_line_strings() process ?file? & ?line? strings

subroutine edit_file_line_strings( edit_line, lower_case_line)

! **********************************************************************

!  edit_file_line_strings() interface

! ----------------------------------------------------------------------

!  the source line to be edited

character( len= *), intent( inout) :: edit_line

!  the source line in lower case to enable searches

character( len= *), intent( inout) :: lower_case_line

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?integer

!  exit: line has any ?name strings replaced with their values

! **********************************************************************

!  edit_file_line_strings() local

! ----------------------------------------------------------------------

!  strings to be edited into the line

   character( len= conversion_len) :: line_number_str

!  find substring

   integer :: search_idx

! **********************************************************************

!  edit_file_line_strings() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  replace ?file? with the current filename

   search_idx = index( lower_case_line, file_str)

   go_file: if( search_idx > 0 )then

      call replace_substring( edit_line, lower_case_line, file_str, trim( current_file% name_str), &
                              search_idx)

      remark_file: if( options% report_extensions )then

         call msg_continue( "edited file: " // trim( current_file% name_str) )

      endif remark_file

   endif go_file

!  replace ?line? with the current line number

   search_idx = index( lower_case_line, line_str)

   go_line: if( search_idx > 0 )then

      write( unit= line_number_str, fmt= conversion_fmt) current_file% lines_transfered

      line_number_str = adjustl( line_number_str)

      call replace_substring( edit_line, lower_case_line, line_str, trim( line_number_str), &
                              search_idx)

      remark_line: if( options% report_extensions )then

         call msg_continue( "edited line: " // trim( line_number_str) )

      endif remark_line

   endif go_line

! ----------------------------------------------------------------------

!  edit_file_line_strings() exit

return

! **********************************************************************

!  edit_file_line_strings()

end subroutine edit_file_line_strings

! **********************************************************************
! **********************************************************************

!  edit_date_time_strings() process ?date? & ?time? strings

subroutine edit_date_time_strings( edit_line, lower_case_line)

! **********************************************************************

!  edit_date_time_strings() interface

! ----------------------------------------------------------------------

!  the source line to be edited

character( len= *), intent( inout) :: edit_line

!  the source line in lower case to enable searches

character( len= *), intent( inout) :: lower_case_line

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?integer

!  exit: line has any ?name strings replaced with their values

! **********************************************************************

!  edit_date_time_strings() local

! ----------------------------------------------------------------------

!  strings to be edited into the line- they are the exact length needed

   character( len= 8) :: today_str

   character( len= 10) :: now_str

!  find substring

   integer :: search_idx

! **********************************************************************

!  edit_date_time_strings() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  replace ?date? with the current date

   search_idx = index( lower_case_line, date_str)

   go_date: if( search_idx > 0 )then

      call date_and_time( date= today_str)

      call replace_substring( edit_line, lower_case_line, date_str, today_str, search_idx)

      remark_date: if( options% report_extensions )then

         call msg_continue( "edited date: " // trim( today_str) )

      endif remark_date

   endif go_date

!  replace ?time? with the current time

   search_idx = index( lower_case_line, time_str)

   go_time: if( search_idx > 0 )then

      call date_and_time( time= now_str)

      call replace_substring( edit_line, lower_case_line, time_str, now_str, search_idx)

      remark_time: if( options% report_extensions )then

         call msg_continue( "edited time: " // trim( now_str) )

      endif remark_time

   endif go_time

! ----------------------------------------------------------------------

!  edit_date_time_strings() exit

return

! **********************************************************************

!  edit_date_time_strings()

end subroutine edit_date_time_strings

! **********************************************************************
! **********************************************************************

!  edit_integer_strings() process ?integer? strings

subroutine edit_integer_strings( edit_line, lower_case_line)

! **********************************************************************

!  edit_integer_strings() interface

! ----------------------------------------------------------------------

!  the source line to be edited

character( len= *), intent( inout) :: edit_line

!  the source line in lower case to enable searches

character( len= *), intent( inout) :: lower_case_line                

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?integer?

!  exit: line has any ?name? strings replaced with their values

! **********************************************************************

!  edit_integer_strings() local

! ----------------------------------------------------------------------

!  string containing integer value

   character( len= conversion_len) :: value_str

!  target string to be replaced

   character( len= target_len) :: search_str

!  length of target string

   integer :: search_len

!  point to integers on symbol list

   type( symbol_t), pointer :: symbol_ptr

!  point to search_str location in line

   integer :: search_idx

! **********************************************************************

!  edit_integer_strings() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  replace ?integer? with the current string value of integer

   nullify( symbol_ptr)

!  test the occurance of each integer on symbol list

   each_integer: do

      call get_next_integer( symbol_ptr)

      end_of_list: if( .not. associated( symbol_ptr) )then

         exit each_integer

      endif end_of_list

!  does ?integer? appear on line

      search_str = arg_key // trim( symbol_ptr% name_str) // arg_key

      search_len = len_trim( search_str)

      search_idx = index( lower_case_line, search_str( : search_len))

!  if found the target, try to replace it with its value

      go_integer: if( search_idx > 0 )then

!  if integer has a value

         defined_integer: if( symbol_ptr% defined )then

            write( unit= value_str, fmt= conversion_fmt) symbol_ptr% integer_value

            value_str = adjustl( value_str)

!  if integer has no value

         else defined_integer

            call msg_quit( "edit integer symbol not defined: " // trim( symbol_ptr% name_str) )

         endif defined_integer

!  go replace the string with its value

         call replace_substring( edit_line, lower_case_line, search_str( : search_len), trim( value_str), &
                                 search_idx)

!  if reporting extensions

         report_integers: if( options% report_extensions )then

            call msg_continue( "edited integer: " // trim( symbol_ptr% name_str) )

         endif report_integers

      endif go_integer

   enddo each_integer

! ----------------------------------------------------------------------

!  edit_integer_strings() exit

return

! **********************************************************************

!  edit_integer_strings()

end subroutine edit_integer_strings

! **********************************************************************
! **********************************************************************

!  edit_logical_strings() process ?logical? strings

subroutine edit_logical_strings( edit_line, lower_case_line)

! **********************************************************************

!  edit_logical_strings() interface

! ----------------------------------------------------------------------

!  the source line to be edited

character( len= *), intent( inout) :: edit_line

!  the source line in lower case to enable searches

character( len= *), intent( inout) :: lower_case_line                     

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?logical

!  exit: line has any ?name strings replaced with their values

! **********************************************************************

!  edit_logical_strings() local

! ----------------------------------------------------------------------

!  string containing logical value

   character( len= conversion_len) :: value_str

!  target string to be replaced

   character( len= target_len) :: search_str

!  length of target string

   integer :: search_len

!  point to logicals on symbol list

   type( symbol_t), pointer :: symbol_ptr

!  point to search_str location in line

   integer :: search_idx

! **********************************************************************

!  edit_logical_strings() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  replace ?logical? with the current string value of logical

   nullify( symbol_ptr)

   each_logical: do

      call get_next_logical( symbol_ptr)

      end_of_list: if( .not. associated( symbol_ptr) )then

         exit each_logical

      endif end_of_list

!  does ?logical? appear on line

      search_str = arg_key // trim( symbol_ptr% name_str) // arg_key

      search_len = len_trim( search_str)

      search_idx = index( lower_case_line, search_str( : search_len))

      go_logical: if( search_idx > 0 )then

!  if logical has a value

         defined_logical: if( symbol_ptr% defined )then

            decode: if( symbol_ptr% logical_value )then

               value_str = true_str

            else decode

               value_str = false_str

            endif decode

            call replace_substring( edit_line, lower_case_line, search_str( : search_len), &
                                    trim( value_str), search_idx)

!  if logical has no value

         else defined_logical

            call msg_quit( "edit logical symbol not defined: " // trim( symbol_ptr% name_str) )

         endif defined_logical

!  if reporting extensions

         report_logicals: if( options% report_extensions )then

            call msg_continue( "edited logical: " // trim( symbol_ptr% name_str) )

         endif report_logicals

      endif go_logical

   enddo each_logical

! ----------------------------------------------------------------------

!  edit_logical_strings() exit

return

! **********************************************************************

!  edit_logical_strings()

end subroutine edit_logical_strings

! **********************************************************************
! **********************************************************************

!  edit_macro_strings() process ?macro? strings

subroutine edit_macro_strings( edit_line, lower_case_line)

! **********************************************************************

!  edit_macro_strings() interface

! ----------------------------------------------------------------------

!  the osurce line to be edited

character( len= *), intent( inout) :: edit_line

!  the source line in lower case to enable searches

character( len= *), intent( inout) :: lower_case_line             

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?macro

!  exit: line has any ?macro strings replaced with their values

! **********************************************************************

!  edit_macro_strings() local

! ----------------------------------------------------------------------

!  string containing macro value

   character( len= target_len) :: search_str

!  length of ?macro?

   integer :: search_len

!  point to ?macro?

   integer :: search_idx

!  argument strings

   character( len= buffer_len) :: actual_args

   character( len= buffer_len) :: value_str

!  scan for macros

   type( symbol_t), pointer :: macro_ptr

!  end of substrings

   integer :: close_paren_idx

   integer :: open_paren_idx

   integer :: value_len

! **********************************************************************

!  edit_macro_strings() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  replace ?macro? with the current string value of macro

   nullify( macro_ptr)

   each_macro: do

      call get_next_macro( macro_ptr)

      end_of_list: if( .not. associated( macro_ptr) )then

         exit each_macro

      endif end_of_list

! ----------------------------------------------------------------------

!  does ?macro? appear on line

      search_str = arg_key // trim( macro_ptr% name_str) // arg_key

      search_len = len_trim( search_str)

      search_idx = index( lower_case_line, search_str( : search_len))

! ----------------------------------------------------------------------

!  if macro appears on line

      found_macro: if( search_idx > 0 )then

!  macro definition has a dummy arg list

         have_arg_list: if( associated( macro_ptr% dummy_args) )then

!  must rebuild the macro value with each new set of actual args

            next_dummy_args: do while( search_idx > 0 )

!  check for actual arg list

               open_paren_idx = search_idx + search_len

               no_actual_args: if( edit_line( open_paren_idx: open_paren_idx) /= open_paren )then

                  call msg_quit( "macro args missing: " // trim( macro_ptr% name_str) )

               endif no_actual_args

!  have an actual arg list, find the close parenthesis

               call seek_close_paren( edit_line, open_paren_idx, close_paren_idx)

               actual_args = lower_case_line( open_paren_idx + 1: close_paren_idx - 1)

!  build the new macro value

               call process_actual_arglist( actual_args( : close_paren_idx - open_paren_idx - 1), &
                                            value_str, macro_ptr)

!  substitute it

               value_len = len_trim( value_str)

!  replace whole "?macro?(args)" with computed macro value

               edit_line = edit_line( : search_idx - 1) // value_str( : value_len) &
                           // edit_line( close_paren_idx + 1: )

               lower_case_line = lower_case_line( : search_idx - 1) // value_str( : value_len) &
                                 // lower_case_line( close_paren_idx + 1: )

!  find the next occurance of ?macro?

               search_idx = index( lower_case_line, trim( search_str( : search_len)) )

            enddo next_dummy_args

! ----------------------------------------------------------------------

!  no arg list so macro value doesn't change

         else have_arg_list

!  insert macro into the line

            value_len = len_trim( macro_ptr% macro_value)

            call replace_substring( edit_line, lower_case_line, search_str( : search_len), &
                                    macro_ptr% macro_value( : value_len), search_idx)

         endif have_arg_list

!  if reporting extensions

         report_macros: if( options% report_extensions )then

            call msg_continue( "edited macro: " // trim( macro_ptr% name_str) )

         endif report_macros

! ----------------------------------------------------------------------

!  done with this macro

      endif found_macro

 ! ----------------------------------------------------------------------

!  go try the next macro

  enddo each_macro

! ----------------------------------------------------------------------

!  edit_macro_strings() exit

return

! **********************************************************************

!  edit_macro_strings()

end subroutine edit_macro_strings

! **********************************************************************
! **********************************************************************

!  process_actual_arglist() process macro dummy arglist strings

subroutine process_actual_arglist( actual_args, value_str, symbol_ptr)

! **********************************************************************

!  process_actual_arglist() interface

! ----------------------------------------------------------------------

!  the comma separated actual args from the macro instance

character( len= *), intent( in) :: actual_args

!  the value of the macro after editing

character( len= *), intent( out) :: value_str

!  a pointer to the macro

type( symbol_t), pointer :: symbol_ptr                              

! **********************************************************************

!  entry: arg_list is an actual argument list
!         macro is a macro variable

!  exit: value_buf has the macro's value with all dummy args replaced by actuals

! **********************************************************************

!  process_actual_arglist() local

! ----------------------------------------------------------------------

!  construct macro text with actual args

   character( len= buffer_len) :: actual_str

   character( len= target_len) :: search_str

   integer :: search_len

! ----------------------------------------------------------------------

!  character pointers

   integer :: this_dummy

   integer :: this_char

   integer :: dummy_idx

   integer :: close_paren_idx

! **********************************************************************

!  process_actual_arglist() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  the value will be edited with the actual args

   value_str = symbol_ptr% macro_value

   actual_str = actual_args

!  put a comma at the end of the actual arglist to process all the actual args

   actual_str( len( actual_args) + 1: len( actual_args) + 1) = comma

!  if checking, verify that all actual args are enclosed in parenthesis

   check_no_parens: if( options% args_in_parens )then

      call verify_actual_parens( actual_str)

   endif check_no_parens

!  loop thru each dummy arg

   each_arg: do this_dummy = 1, size( symbol_ptr% dummy_args)

!  prepare the ?dummy? string

      search_str = arg_key // trim( symbol_ptr% dummy_args( this_dummy) ) // arg_key

      search_len = len_trim( search_str)

!  scan through the actual arguments string to find a comma outside parenthesis

      this_char = 1

      each_char: do

!  find a comma outside parenthesis

         find_actual: select case( actual_str( this_char: this_char))

!  at open paren, skip to matching paren

         case( open_paren) find_actual

            call seek_close_paren( actual_str, this_char, close_paren_idx)

            this_char = close_paren_idx + 1

!  actual argument is isolated before comma outside parenthesis

         case( comma) find_actual

            dummy_idx = index( value_str, search_str( : search_len) )

            substitute: if( dummy_idx > 0 )then

               call replace_substring( lower_case_str= value_str, search_str= search_str( : search_len), &
                                       replace_str= actual_str( : this_char - 1), first_idx= dummy_idx)

               dummy_idx = index( value_str, search_str( : search_len))

            endif substitute

            actual_str = actual_str( this_char + 1: )

            exit each_char

!  otherwise, keep checking characters

         case default find_actual

            this_char = this_char + 1

         end select find_actual

      enddo each_char

   enddo each_arg

! ----------------------------------------------------------------------

!  process_actual_arglist() exit

return

! **********************************************************************

!  process_actual_arglist()

end subroutine process_actual_arglist

! **********************************************************************
! **********************************************************************

!  verify_actual_parens() verify that all actual arguments have enclosing parenthesis

subroutine verify_actual_parens( actual_args)

! **********************************************************************

!  verify_actual_parens() interface

! ----------------------------------------------------------------------

!  the comma separated actual args from the macro instance

character( len= *), intent( in) :: actual_args

! **********************************************************************

!  entry: actual_args

!  exit: if all actual args are enclosed in paren, or need not be

! **********************************************************************

!  verify_actual_parens() local

! ----------------------------------------------------------------------

!  construct macro text with actual args

   character( len= buffer_len) :: actual_str

   character( len= buffer_len) :: one_actual_arg

! ----------------------------------------------------------------------

!  character pointers

   integer :: this_char

   logical :: all_chars_ok

   integer :: end_skip_idx

! **********************************************************************

!  verify_actual_parens() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  the value will be edited with the actual args

   actual_str = actual_args

!  loop thru each dummy arg

   each_arg: do while( actual_str /= blank)

!  scan through the actual arguments string to find a comma outside parenthesis

      all_chars_ok = .true.

      this_char = 1

      each_char: do

!  find a comma outside parenthesis

         check_char: select case( actual_str( this_char: this_char))

!  at open paren, skip to matching paren

         case( open_paren) check_char

            call seek_close_paren( actual_str, this_char, end_skip_idx)

            this_char = end_skip_idx + 1

!  at quote, skip to matching quote

         case( single_quote, double_quote) check_char

            call seek_match_quote( actual_str, this_char, end_skip_idx)

            this_char = end_skip_idx

!  at allowable character, skip it

         case( 'a': 'z', '0': '9', '_', '%', blank) check_char

            this_char = this_char + 1

!  actual argument is isolated at comma outside parenthesis

         case( comma) check_char

            one_actual_arg = actual_str( : this_char - 1)

            actual_str = actual_str( this_char + 1: )

            exit each_char

!  otherwise, keep checking characters

         case default check_char

            all_chars_ok = .false.

            this_char = this_char + 1

         end select check_char

      enddo each_char

      got_bad_char: if( .not. all_chars_ok )then

         call msg_continue( 'actual argument may need parenthesis: ' // trim( one_actual_arg))

      endif got_bad_char

   enddo each_arg

! ----------------------------------------------------------------------

!  verify_actual_parens() exit

return

! **********************************************************************

!  verify_actual_parens()

end subroutine verify_actual_parens

! **********************************************************************
! **********************************************************************

!  wrap_source_line() ensure lines are not too long

subroutine wrap_source_line( wrap_line)

! **********************************************************************

!  wrap_source_line() interface

! ----------------------------------------------------------------------

!  the line to be wrapped

character( len= *), target, intent( inout) :: wrap_line           

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) more than 132 characters

!  exit: line has continuations written and fewer than 132 source lines

! **********************************************************************

!  wrap_source_line() constants

! ----------------------------------------------------------------------

!  break after one of these

   character( len= *), parameter :: operators = '+-*/=,()%'

   character( len= *), parameter :: break_point = operators // white_space

   character( len= *), parameter :: blank_str = repeat( blank, fixed_format_len)

!  start in the right column for fixed format

   integer, parameter :: start_col = 5

! **********************************************************************

!  wrap_source_line() local

! ----------------------------------------------------------------------

!  length of source line

   integer :: output_len

!  location to wrap lines

   integer :: break_idx

!  flag wrapping

   logical :: did_wrap

! **********************************************************************

!  wrap_source_line() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  initialize

   did_wrap = .false.

   output_len = len_trim( wrap_line)

!  while line is too long

   wrap_lines: do while( output_len > options% wrap_length )

!  break at a break point before the maximum length

      break_idx = scan( wrap_line( : options% wrap_length - 1), break_point, back= .true.)

      no_break: if( break_idx == 0 )then

         break_idx = options% wrap_length

      endif no_break

!  process fixed format differently

      fix_length: if( options% wrap_length == fixed_format_len )then ! fixed or free format wrapping

!  fixed format line up to the breakpoint, then blanks to column 72, then & in column 73

         line = wrap_line( : break_idx - 1) // blank_str( break_idx: fixed_format_len) // continuation

!  blanks up to column 6, then &, then the rest of the line

         wrap_line = blank_str( : start_col) // continuation // line( break_idx: )

!  process free format

      else fix_length

!  free format line up to the breakpoint, then the continuation character

         line = wrap_line( : break_idx - 1) // continuation

! the &, then the rest of the line

         wrap_line = continuation // line( break_idx: )

      endif fix_length

!  write lines as they are made

      call write_source_line( output_file)

      output_len = len_trim( line)

      did_wrap = .true.

!  while line is too long

   enddo wrap_lines

! ----------------------------------------------------------------------

!  if reporting extensions

   extensions: if( did_wrap .and. options% report_extensions )then

      call msg_continue( "wrapped source line")

   endif extensions

! ----------------------------------------------------------------------

!  wrap_source_line() exit

return

! **********************************************************************

!  wrap_source_line()

end subroutine wrap_source_line

! **********************************************************************
! **********************************************************************

!  blank_compress_lower_case() blank compress or convert to lower case

subroutine blank_compress_lower_case( out_str, in_str)

! **********************************************************************

!  blank_compress_lower_case() interface

! ----------------------------------------------------------------------

!  a coco line with blanks, quoted strings and comments

character( len= *), intent( in) :: in_str

!  a blank compressed lower case coco statement

character( len= *), intent( out) :: out_str                         

! **********************************************************************

!  entry: in_str is a coco line

!  exit: out_str is a coco statement, possibly complete only up to hwm

! **********************************************************************

!  blank_compress_lower_case() local

! ----------------------------------------------------------------------

!  quote used for the current quoted string

   character( len= 1), save :: quote

!  length of in_str

   integer :: in_str_len

!  input pointer reset for each line

   integer :: in_idx

!  output pointer reset for each statement

   integer, save :: out_idx = 0

!  character pointer reset by each intrinsic use

   integer :: char_idx

! **********************************************************************

!  blank_compress_lower_case() text

continue

! ----------------------------------------------------------------------

!  null input signals end-of-file reached on input file

   cleanup: if( in_str == null_string )then

      out_idx = 0

      return

   endif cleanup

! ----------------------------------------------------------------------

!  initialize line length

   in_str_len = len_trim( in_str)

!  setup pointers

   initialize: if( out_idx == 0 )then

      in_idx = 1

      out_idx = 1

      out_str = null_string

      quote = null_string

   else initialize

      char_idx = verify( in_str, white_space)

!  check whether first character is a continuation

      skip_contin: if( in_str( char_idx: char_idx) == continuation )then

         in_idx = char_idx + 1

      else skip_contin

         in_idx = char_idx

      endif skip_contin

   endif initialize

! **********************************************************************

!  scan each character until end of input string

   scan_line: do while( in_idx <= in_str_len)

! ----------------------------------------------------------------------

!  if in quoted string

      char_literal: select case( quote)

! ----------------------------------------------------------------------

      case( single_quote, double_quote) char_literal

!  if found matching quote

         end_of_string: if( in_str( in_idx: in_idx) == quote )then

!  out of string so set quote to null

            quote = null_string

         endif end_of_string

!  copy the character

         out_str( out_idx: out_idx) = in_str( in_idx: in_idx)

!  update the pointers

         in_idx = in_idx + 1

         out_idx = out_idx + 1

! ----------------------------------------------------------------------

!  not in quoted string

      case default char_literal

!  white space is not copied

         skip_ws: select case( in_str( in_idx: in_idx) )

!  blanks or tabs

         case( blank, tab) skip_ws

            in_idx = in_idx + 1

!  all others

         case default skip_ws

!  check for special characters
 
            spec_char: select case( in_str( in_idx: in_idx) )

!  found quoted string

            case( single_quote, double_quote) spec_char

               quote = in_str( in_idx: in_idx)

!  found coco comment

            case( comment) spec_char

               exit scan_line

            end select spec_char

! ----------------------------------------------------------------------

!  copy non-blank characters

            copy_char: select case( in_str( in_idx: in_idx) )

!  copy the upper case character as lower case characters

            case( 'A': 'Z') copy_char

               out_str( out_idx: out_idx) = achar( iachar( in_str( in_idx: in_idx)) + change_case)

            case default copy_char

!  copy the other non-blank character unchanged

               out_str( out_idx: out_idx) = in_str( in_idx: in_idx)

            end select copy_char

!  update pointers

            in_idx = in_idx + 1

            out_idx = out_idx + 1

         end select skip_ws

! ----------------------------------------------------------------------

      end select char_literal

!  process next character

   enddo scan_line

! ----------------------------------------------------------------------

!  check whether last character is continuation

   line_complete: if( out_str( out_idx - 1: out_idx - 1) == continuation )then

!  next line is a continuation line

      out_idx = out_idx - 1

   else line_complete

!  next line is an initial line

      out_idx = 0

   endif line_complete

! ----------------------------------------------------------------------

!  blank_compress_lower_case() exit

return

! **********************************************************************

!  blank_compress_lower_case()

end subroutine blank_compress_lower_case

! **********************************************************************
! **********************************************************************

!  process_coco_statement() process a coco directive

subroutine process_coco_statement( coco_stmt)

! **********************************************************************

!  process_coco_statement() interface

! ----------------------------------------------------------------------

!  the coco statement to be processed

character( len= *), intent( in) :: coco_stmt

! **********************************************************************

!  entry: coco_stmt is a blank_compress_lower_case coco directive past the coco key
!         "stop..." | "message..." | "if..." | "elseif..." | "else..." |
!         "endif..." | "integer..." | "logical..." | "assert..." | "<name>=..." |
!         "text..." | "copy..." | "dump" | "options" | "report" |
!         "assertif(" | "copyif(" | "ifndef(" | "ifdef(" | "undefine:" | "doc" |
!         "output"

!  exit: the directive is processed or error exit

!  If a directive might have something after the keyword, the keyword
!  match is checked by "keyword( : len( keyword) ) == string", otheriwse,
!  if the directive must not have anything after the keyword, the
!  keyword match is checked by "keyword == string".  Thus, a directive
!  with unknown nonblank characters after the keyword is an unknown directive.

! **********************************************************************

!  process_coco_statement() local

! ----------------------------------------------------------------------

!  point to location of symbol

   type( symbol_t), pointer :: symbol_ptr

!  possible index of equals

   integer :: eq_idx

!  expression string is after the equals

   integer :: expr_idx

! **********************************************************************

!  process_coco_statement() text

continue

! ----------------------------------------------------------------------

!  detect assignment statements assigning to variables named by keywords

   nullify( symbol_ptr)

   eq_idx = scan( coco_stmt( : symbol_name_len + len( equals)), equals)

   got_equals: if( eq_idx > 0 )then

      call seek_symbol_name( coco_stmt( : eq_idx - 1), symbol_ptr)

   endif got_equals

! ----------------------------------------------------------------------

!  which directive?

! ----------------------------------------------------------------------

!  assignment directive

   which_directive: if( associated( symbol_ptr) )then

!  up to the equals must be a declared name

      expr_idx = eq_idx + len( equals)

!  must be an integer or logical variable

      integer_or_logical: select case( symbol_ptr% type_code)
         
      case( type_integer) integer_or_logical

         call process_integer_assignment( coco_stmt( expr_idx: ), symbol_ptr)

      case( type_logical) integer_or_logical

         call process_logical_assignment( coco_stmt( expr_idx: ), symbol_ptr)

      case default integer_or_logical

         call msg_quit( "variable must be an integer or a logical: " // trim( symbol_ptr% name_str) )
         
      end select integer_or_logical

      nullify( symbol_ptr)

! ----------------------------------------------------------------------

!  stop directive

   elseif( coco_stmt == stop_str )then which_directive

      call process_stop_directive( coco_stmt( len( stop_str) + 1: ) )

! ----------------------------------------------------------------------

!  message directive

   elseif( coco_stmt( : len( message_str)) == message_str )then which_directive

      call process_message_directive( coco_stmt( len( message_str) + 1: ) )

! ----------------------------------------------------------------------

!  if directive

   elseif( coco_stmt( : len( if_str)) == if_str )then which_directive

      call process_if_directive( coco_stmt( len( if_str) + 1: ) )

! ----------------------------------------------------------------------

!  elseif directive

   elseif( coco_stmt( : len( elseif_str)) == elseif_str )then which_directive

      call process_elseif_directive( coco_stmt( len( elseif_str) + 1: ) )

! ----------------------------------------------------------------------

!  else directive

   elseif( coco_stmt( : len( else_str)) == else_str )then which_directive

      call process_else_directive( coco_stmt( len( else_str) + 1: ) )

! ----------------------------------------------------------------------

!  endif directive

   elseif( coco_stmt( : len( endif_str)) == endif_str )then which_directive

      call process_endif_directive( coco_stmt( len( endif_str) + 1: ) )

! ----------------------------------------------------------------------

!  ifdef directive

   elseif( coco_stmt( : len( ifdef_str)) == ifdef_str )then which_directive

      call process_ifdef_directive( coco_stmt( len( ifdef_str) + 1: ) )

! ----------------------------------------------------------------------

!  ifndef directive

   elseif( coco_stmt( : len( ifndef_str)) == ifndef_str )then which_directive

      call process_ifndef_directive( coco_stmt( len( ifndef_str) + 1: ) )

! ----------------------------------------------------------------------

!  integer declaration

   elseif( coco_stmt( : len( integer_str)) == integer_str )then which_directive

      call process_integer_declaration( coco_stmt( len( integer_str) + 1: ))

! ----------------------------------------------------------------------

!  integer constant declaration

   elseif( coco_stmt( : len( integer_constant_str)) == integer_constant_str )then which_directive

      call process_integer_constant( coco_stmt( len( integer_constant_str) + 1: ))

! ----------------------------------------------------------------------

!  logical declaration

   elseif( coco_stmt( : len( logical_str)) == logical_str )then which_directive

      call process_logical_declaration( coco_stmt( len( logical_str) + 1: ))

! ----------------------------------------------------------------------

!  logical constant declaration

   elseif( coco_stmt( : len( logical_constant_str)) == logical_constant_str )then which_directive

      call process_logical_constant( coco_stmt( len( logical_constant_str) + 1: ))

! ----------------------------------------------------------------------

!  macro declaration

   elseif( coco_stmt( : len( macro_str)) == macro_str )then which_directive

      call process_macro_declaration( coco_stmt( len( macro_str) + 1: ))

! ----------------------------------------------------------------------

!  undefine declaration

   elseif( coco_stmt( : len( undefine_str)) == undefine_str )then which_directive

      call process_undefine_directive( coco_stmt( len( undefine_str) + 1: ))

! ----------------------------------------------------------------------

!  assertif directive -- assertif must come before assert due to leading substring

   elseif( coco_stmt( : len( assertif_str)) == assertif_str )then which_directive

      call process_assertif_directive( coco_stmt( len( assertif_str): ))

! ----------------------------------------------------------------------

!  assert directive

   elseif( coco_stmt( : len( assert_str)) == assert_str )then which_directive

      call process_assert_directive( coco_stmt( len( assert_str) + 1: ))

! ----------------------------------------------------------------------

!  dump directive

   elseif( coco_stmt == dump_str )then which_directive

      call process_dump_directive

! ----------------------------------------------------------------------

!  options directive

   elseif( coco_stmt == options_str )then which_directive

      call write_options

! ----------------------------------------------------------------------

!  report directive

   elseif( coco_stmt == report_str )then which_directive

      call write_report

! ----------------------------------------------------------------------

!  text directive

   elseif( coco_stmt( : len( text_str)) == text_str )then which_directive

      call process_text_directive( coco_stmt( len( text_str) + 1: ))

! ----------------------------------------------------------------------

!  copy directive

   elseif( coco_stmt( : len( copy_str)) == copy_str )then which_directive

      call process_copy_directive( coco_stmt( len( copy_str) + 1: ))

! ----------------------------------------------------------------------

!  copyif directive

   elseif( coco_stmt( : len( copyif_str)) == copyif_str )then which_directive

      call process_copyif_directive( coco_stmt( len( copyif_str): ))

! ----------------------------------------------------------------------

!  doc directive

   elseif( coco_stmt( : len( doc_str)) == doc_str )then which_directive

      call process_doc_directive( coco_stmt( len( doc_str) + 1: ))

! ----------------------------------------------------------------------

!  output directive

   elseif( coco_stmt( : len( output_str)) == output_str )then which_directive

      call process_output_directive( coco_stmt( len( output_str) + 1: ))

! ----------------------------------------------------------------------

!  cannot process this directive

   else which_directive

      call msg_quit( "unknown coco directive: " // trim( coco_stmt))

! ----------------------------------------------------------------------

!  which directive?

   endif which_directive

! ----------------------------------------------------------------------

!  process_coco_statement() exit

return

! **********************************************************************

!  process_coco_statement()

end subroutine process_coco_statement

! **********************************************************************
! **********************************************************************

!  process_integer_assignment() process a coco stop directive

subroutine process_integer_assignment( assign_dir, integer_ptr)

! **********************************************************************

!  process_integer_assignment() interface

! ----------------------------------------------------------------------

!  the integer assignment directive

character( len= *), intent( in) :: assign_dir

!  a pointer to the integer symbol

type( symbol_t), pointer :: integer_ptr

! **********************************************************************

!  entry: stop_dir is blank_compress_lower_case coco stop directive, past the coco key word

!  exit: coco processing stops

! **********************************************************************

!  process_integer_assignment() local

! ----------------------------------------------------------------------

!  value is an intent( out) variable

   integer :: value

! **********************************************************************

!  process_integer_assignment() text

continue

! ----------------------------------------------------------------------

!  process assignment directive if on an active line

   active_line: if( if_construct% now_selected ) then

!  do not allow redefinition of constants

      redefine_constant: if( integer_ptr% constant )then

         call msg_quit( "attempt to redefine a constant: " // trim( integer_ptr% name_str) )

      endif redefine_constant

!  assign the value

      call eval_int_expr( assign_dir, value)

      integer_ptr% integer_value = value

      integer_ptr% defined = .true.

   endif active_line

! ----------------------------------------------------------------------

!  process_integer_assignment() exit

return

! **********************************************************************

!  process_integer_assignment()

end subroutine process_integer_assignment

! **********************************************************************
! **********************************************************************

!  process_logical_assignment() process a coco stop directive

subroutine process_logical_assignment( assign_dir, logical_ptr)

! **********************************************************************

!  process_logical_assignment() interface

! ----------------------------------------------------------------------

!  the logical assignment directive

character( len= *), intent( in) :: assign_dir

!  a pointer to the logical symbol

type( symbol_t), pointer :: logical_ptr

! **********************************************************************

!  entry: stop_dir is blank_compress_lower_case coco stop directive, past the coco key word

!  exit: coco processing stops

! **********************************************************************

!  process_logical_assignment() local

! ----------------------------------------------------------------------

!  value is an intent( out) variable

   logical :: value

! **********************************************************************

!  process_logical_assignment() text

continue

! ----------------------------------------------------------------------

!  process stop directive if on an active line

   active_line: if( if_construct% now_selected )then

!  do not allow redefinition of constants

      redefine_constant: if( logical_ptr% constant )then

         call msg_quit( "attempt to redefine a constant: " // trim( logical_ptr% name_str) )

      endif redefine_constant

!  assign the value

      call eval_log_expr( assign_dir, value)

      logical_ptr% logical_value = value

      logical_ptr% defined = .true.

   endif active_line

! ----------------------------------------------------------------------

!  process_logical_assignment() exit

return

! **********************************************************************

!  process_logical_assignment()

end subroutine process_logical_assignment

! **********************************************************************
! **********************************************************************

!  process_stop_directive() process a coco stop directive

subroutine process_stop_directive( stop_dir)

! **********************************************************************

!  process_stop_directive() interface

! ----------------------------------------------------------------------

!  the stop directive

character( len= *), intent( in) :: stop_dir                         

! **********************************************************************

!  entry: stop_dir is blank_compress_lower_case coco stop directive, past the coco key word

!  exit: coco processing stops

! **********************************************************************

!  process_stop_directive() text

continue

! ----------------------------------------------------------------------

!  process stop directive if on an active line

   active_line: if( if_construct% now_selected )then

      verbose_output: if( options% verbose_mode )then

         call msg_continue( "coco stop directive encountered: " // trim( stop_dir))

      endif verbose_output

      output_file% line = stop_dir

      call write_coco_line( output_file)

      stop 'coco stop directive encountered'

   endif active_line

! ----------------------------------------------------------------------

!  process_stop_directive() exit

return

! **********************************************************************

!  process_stop_directive()

end subroutine process_stop_directive

! **********************************************************************
! **********************************************************************

!  process_message_directive() process a coco message directive

subroutine process_message_directive( message_dir)

! **********************************************************************

!  process_message_directive() interface

! ----------------------------------------------------------------------

!  the message directive

character( len= *), intent( in) :: message_dir                       

! **********************************************************************

!  entry: message_dir is blank_compress_lower_case coco message directive, past the message key word

!  exit: message is written to error unit

! **********************************************************************

!  process_message_directive() local

! ----------------------------------------------------------------------

   character( len= buffer_len) :: msg_buffer

   integer :: in_idx

   integer :: out_idx

   integer :: comma_idx

   integer :: quoted_len

   integer :: unquoted_len

   integer :: integer_value

   integer :: istat

   logical :: logical_value

   logical :: is_integer

   character( len= buffer_len) :: expr_str

   character( len= conversion_len) :: conversion_str

! **********************************************************************

!  process_message_directive() text

continue

! ----------------------------------------------------------------------

!  process if on active line

   active_line: if( if_construct% now_selected )then

! ----------------------------------------------------------------------

!  initialize

      in_idx = 1
      out_idx = 1

      msg_buffer = blank

!  loop thru message list items

      list_items: do while( in_idx <= len_trim( message_dir) )

! ----------------------------------------------------------------------

!  a list item can be a quoted string or an expression

         string_expr: select case( message_dir( in_idx: in_idx) )

! ----------------------------------------------------------------------

!  process quoted strings

         case( single_quote, double_quote) string_expr

!  try to unquote the string

            call unquote_string( message_dir( in_idx: ), msg_buffer( out_idx: ), quoted_len, unquoted_len)

!  if the matching quote is found within the string

            got_string: if( quoted_len <= len( message_dir( in_idx: )) )then

!  found quote, update the character pointers

               in_idx = in_idx + quoted_len + index( message_dir( in_idx + quoted_len - 1: ), comma) - 1

               out_idx = out_idx + unquoted_len - 1

!  found not quote, complain and quit

            else got_string

               call msg_quit( "bad message string: " // trim( message_dir) )

            endif got_string

! ----------------------------------------------------------------------

!  process expressions

         case( 'a': 'z', '0':'9', dot, plus, minus, open_paren) string_expr

!  expression ends at a comma

            comma_idx = scan( message_dir( in_idx: ), comma)

!  find the comma or the end of the expression

            end_of_string: if( comma_idx == 0 )then

               comma_idx = len_trim( message_dir) + 1

            else end_of_string

               comma_idx = in_idx + comma_idx - 1

            endif end_of_string

!  encode integer or logical

            call integer_or_logical( message_dir( in_idx: comma_idx - 1), is_integer)

!  an integer expression

            int_log: if( is_integer )then

               expr_str = message_dir( in_idx: comma_idx - 1)

               call eval_int_expr( expr_str, integer_value)

               write( unit= conversion_str, fmt= conversion_fmt, iostat= istat) integer_value

!  trap internal write errors

               encode: if( istat > 0 )then

                  call msg_quit( "can't encode: " // message_dir( in_idx: comma_idx - 1) )

               endif encode

               msg_buffer( out_idx: ) = adjustl( conversion_str)

!  a logical expression

            else int_log

               expr_str = message_dir( in_idx: comma_idx - 1)

               call eval_log_expr( expr_str, logical_value)

               t_or_f: if( logical_value )then

                  msg_buffer( out_idx: ) = '.true.'

               else t_or_f

                  msg_buffer( out_idx: ) = '.false.'

               endif t_or_f

            endif int_log

!  update pointers and add to output buffer

            adjust: if( msg_buffer( out_idx: out_idx) == blank )then ! space out output

               msg_buffer( out_idx + 1: ) = adjustl( msg_buffer( out_idx + 1: ) )

            endif adjust

            in_idx = comma_idx + 1

            out_idx = len_trim( msg_buffer) + 2

! ----------------------------------------------------------------------

!  list item isn't a string, a symbol or a literal

         case default string_expr

            call msg_quit( "bad message list item: " // message_dir( in_idx: ) )

! ----------------------------------------------------------------------

         end select string_expr

!  loop thru message list items

      enddo list_items

! ----------------------------------------------------------------------

!  make the message available

      verbose_output: if( options% verbose_mode )then

         call msg_continue( msg_buffer( : out_idx) )

      endif verbose_output

      call wrap_source_line( msg_buffer( : out_idx) )

      line = msg_buffer

      call write_coco_line( output_file)

   endif active_line

! ----------------------------------------------------------------------

!  process_message_directive() exit

return

! **********************************************************************

!  process_message_directive()

end subroutine process_message_directive

! **********************************************************************
! **********************************************************************

!  %%% process declarations of integers, logicals, macros and texts

! **********************************************************************
! **********************************************************************

!  get_symbol_name() verify symbol name and determine its length

subroutine get_symbol_name( decl_stmt, symbol_name, name_len)

! **********************************************************************

!  get_symbol_name() interface

! ----------------------------------------------------------------------

!  a declaration statement with a symbol name

character( len= *), intent( in) :: decl_stmt

!  the name of the symbol

character( len= *), intent( out) :: symbol_name

!  the length of the symbol name

integer, intent( out) :: name_len

! **********************************************************************

!  entry: decl_stmt is blank_compress_lower_case declaration statement past the double colon
!         "name" | "name=..."

!  exit: a valid symbol name and its length or error exit

! **********************************************************************

!  get_symbol_name() constants

! ----------------------------------------------------------------------

!  characters which must end a symbol name

   character( len= *), parameter :: end_of_name = equals // blank

! **********************************************************************

!  get_symbol_name() local

! ----------------------------------------------------------------------

!  pointers to characters in decl_stmt

   integer :: char_idx

! **********************************************************************

!  get_symbol_name() text

continue

! ----------------------------------------------------------------------

!  look for equals following separator

   char_idx = scan( decl_stmt, end_of_name)

   name_error: if( char_idx == 0 )then

      call msg_quit( "can't find name in declaration: " // trim( decl_stmt))

   endif name_error

!  length of name is one less than first character past name

   name_len = char_idx - 1

! ----------------------------------------------------------------------

!  check that initial character is alphabetic

   char_idx = verify( decl_stmt( 1: 1), alpha_chars)

   initial_ok: if( char_idx /= 0 )then

      call msg_quit( "illegal initial character in symbol name: " // trim( decl_stmt) )

   endif initial_ok

!  check that following characters are legal

   char_idx = verify( decl_stmt( 2: name_len), alphanum_chars)

   name_ok: if( char_idx /= 0 )then

      call msg_quit( "illegal character in symbol name: " // trim( decl_stmt) )

   endif name_ok

!  return name

   symbol_name = decl_stmt( : name_len)

! ----------------------------------------------------------------------

!  get_symbol_name() exit

return

! **********************************************************************

!  get_symbol_name()

end subroutine get_symbol_name

! **********************************************************************
! **********************************************************************

!  get_macro_name() verify macro name and determine its length

subroutine get_macro_name( decl_stmt, macro_name, name_len)

! **********************************************************************

!  get_macro_name() interface

! ----------------------------------------------------------------------

!  the directive containing the macro

character( len= *), intent( in) :: decl_stmt

!  the name of the macro

character( len= *), intent( out) :: macro_name

!  the length of the macro name

integer, intent( out) :: name_len

! **********************************************************************

!  entry: decl_stmt is blank_compress_lower_case declaration statement past the double colon
!         "name=..." | "name(..."

!  exit: name is valid and its length is known or error exit

! **********************************************************************

!  get_macro_name() constants

! ----------------------------------------------------------------------

!  equals or open parenthesis may end a macro name

   character( len= *), parameter :: end_of_name = equals // open_paren

! **********************************************************************

!  get_macro_name() local

! ----------------------------------------------------------------------

!  pointers to characters in decl_stmt

   integer :: char_idx

! **********************************************************************

!  get_macro_name() text

continue

! ----------------------------------------------------------------------

!  look for equals following separator

   char_idx = scan( decl_stmt, end_of_name)

!  if no equals or open paren found

   no_eq_op: if( char_idx == 0 )then

     call msg_quit( "can't find name in macro declaration: " // trim( decl_stmt))

   endif no_eq_op

   name_len = char_idx - 1

! ----------------------------------------------------------------------

!  check that initial character is alphabetic

   char_idx = verify( decl_stmt( 1: 1), alpha_chars)

   initial_ok: if( char_idx /= 0 )then

      call msg_quit( "illegal initial character in macro name: " // trim( decl_stmt) )

   endif initial_ok

!  check that following characters are legal

   char_idx = verify( decl_stmt( 2: name_len), alphanum_chars)

   name_ok: if( char_idx /= 0 )then

      call msg_quit( "illegal character in macro name: " // trim( decl_stmt) )

   endif name_ok

   macro_name = decl_stmt( : name_len)

! ----------------------------------------------------------------------

!  get_macro_name() exit

return

! **********************************************************************

!  get_macro_name()

end subroutine get_macro_name

! **********************************************************************
! **********************************************************************

!  get_text_name() verify text name and determine its length

subroutine get_text_name( decl_stmt, text_name, name_len)

! **********************************************************************

!  get_text_name() interface

! ----------------------------------------------------------------------

!  the statement containing the text name

character( len= *), intent( in) :: decl_stmt

!  the text name

character( len= *), intent( out) :: text_name

!  the length of the text name

integer, intent( out) :: name_len                                   

! **********************************************************************

!  entry: decl_stmt is blank_compress_lower_case declaration statement past the double colon
!         "name" | "name(..."

!  exit: name is valid and its length is known or error exit

! **********************************************************************

!  get_text_name() constants

! ----------------------------------------------------------------------

!  blank or open parenthesis may end a name

   character( len= *), parameter :: end_of_name = blank // open_paren

! **********************************************************************

!  get_text_name() local

! ----------------------------------------------------------------------

!  pointers to characters in decl_stmt

   integer :: char_idx

! **********************************************************************

!  get_text_name() text

continue

! ----------------------------------------------------------------------

!  look for equals following separator

   char_idx = scan( decl_stmt, end_of_name)

!  if no equals found

   no_eq_op: if( char_idx == 0 )then

     call msg_quit( "can't find name in macro statement: " // trim( decl_stmt))

   endif no_eq_op

   name_len = char_idx - 1

! ----------------------------------------------------------------------

!  check that initial character is alphabetic

   char_idx = verify( decl_stmt( 1: 1), alpha_chars)

   initial_ok: if( char_idx /= 0 )then

      call msg_quit( "illegal initial character in text name: " // trim( decl_stmt) )

   endif initial_ok

!  check that following characters are legal

   char_idx = verify( decl_stmt( 2: name_len), alphanum_chars)

   name_ok: if( char_idx /= 0 )then

      call msg_quit( "illegal character in text name: " // trim( decl_stmt) )

   endif name_ok

   text_name = decl_stmt( : name_len)

! ----------------------------------------------------------------------

!  get_text_name() exit

return

! **********************************************************************

!  get_text_name()

end subroutine get_text_name

! **********************************************************************
! **********************************************************************

!  process_integer_declaration() process integer declarations

subroutine process_integer_declaration( integer_stmt)

! **********************************************************************

!  process_integer_declaration() interface

! ----------------------------------------------------------------------

!  the statement containing the declaration

character( len= *), intent( in) :: integer_stmt                    

! **********************************************************************

!  entry: int_stmt is blank_compress_lower_case integer declaration past the integer keyword
!         "::..." | ",parameter::..."

!  exit: int_stmt is processed or error exit

! **********************************************************************

!  process_integer_declaration() constants

! ----------------------------------------------------------------------

!  mark the end of a definition

   character( len= *), parameter :: end_of_def = comma // blank

! **********************************************************************

!  process_integer_declaration() local

! ----------------------------------------------------------------------

!  string containing a single symbol declaration symbol

   character( len= buffer_len) :: symbol_str

!  name of symbol

   character( len= symbol_name_len) :: symbol_name

!  results of decoding statement

   integer :: symbol_len

! ----------------------------------------------------------------------

!  point to next character to be decoded

   integer :: next_char

   integer :: def_len

! **********************************************************************

!  process_integer_declaration() text

continue

   next_char = 1

! ----------------------------------------------------------------------

!  if active line, process the declaration

   active_line: if( if_construct% now_selected )then

!  extract all symbols on directive

      all_symbols: do

!  one symbol at a time to the symbol string

         def_len = scan( integer_stmt( next_char: ), end_of_def) + next_char - 1

         symbol_str = integer_stmt( next_char: def_len - 1)

!  extract symbol name

         call get_symbol_name( symbol_str, symbol_name, symbol_len)

!  store symbol in symbol list

         call add_integer( symbol_str( symbol_len + 1: def_len - 1), symbol_name, .false.)

!  comma separates symbols, blank is end of statement

         all_done: if( integer_stmt( def_len: def_len) == blank )then

            exit all_symbols

         endif all_done

!  move to next symbol

         next_char = def_len + 1

!  extract all symbols on directive

      enddo all_symbols

!  if active line, process the declaration

   endif active_line

! ----------------------------------------------------------------------

!  process_integer_declaration() exit

return

! **********************************************************************

!  process_integer_declaration()

end subroutine process_integer_declaration

! **********************************************************************
! **********************************************************************

!  process_integer_constant() process integer declarations

subroutine process_integer_constant( integer_stmt)

! **********************************************************************

!  process_integer_constant() interface

! ----------------------------------------------------------------------

!  the statement containing the declaration

character( len= *), intent( in) :: integer_stmt                    

! **********************************************************************

!  entry: int_stmt is blank_compress_lower_case integer declaration past the integer keyword
!         "::..." | ",parameter::..."

!  exit: int_stmt is processed or error exit

! **********************************************************************

!  process_integer_constant() constants

! ----------------------------------------------------------------------

!  mark the end of a definition

   character( len= *), parameter :: end_of_def = comma // blank

! **********************************************************************

!  process_integer_constant() local

! ----------------------------------------------------------------------

!  string containing a single symbol declaration symbol

   character( len= buffer_len) :: symbol_str

!  name of symbol

   character( len= symbol_name_len) :: symbol_name

!  results of decoding statement

   integer :: symbol_len

! ----------------------------------------------------------------------

!  point to next character to be decoded

   integer :: next_char

   integer :: def_len

! **********************************************************************

!  process_integer_constant() text

continue

   next_char = 1

! ----------------------------------------------------------------------

!  if active line, process the declaration

   active_line: if( if_construct% now_selected )then

!  extract all symbols on directive

      all_symbols: do

!  one symbol at a time to the symbol string

         def_len = scan( integer_stmt( next_char: ), end_of_def) + next_char - 1

         symbol_str = integer_stmt( next_char: def_len - 1)

!  extract symbol name

         call get_symbol_name( symbol_str, symbol_name, symbol_len)

!  store symbol in symbol list

         call add_integer( symbol_str( symbol_len + 1: def_len - 1), symbol_name, .true.)

!  comma separates symbols, blank is end of statement

         all_done: if( integer_stmt( def_len: def_len) == blank )then

            exit all_symbols

         endif all_done

!  move to next symbol

         next_char = def_len + 1

!  extract all symbols on directive

      enddo all_symbols

!  if active line, process the declaration

   endif active_line

! ----------------------------------------------------------------------

!  process_integer_constant() exit

return

! **********************************************************************

!  process_integer_constant()

end subroutine process_integer_constant

! **********************************************************************
! **********************************************************************

!  add_integer() store integer declaration in symbol table

subroutine add_integer( int_decl_str, symbol_name, is_const)

! **********************************************************************

!  add_integer() interface

! ----------------------------------------------------------------------

!  the statement containing the declaration

character( len= *), intent( in) :: int_decl_str

!  the symbol name

character( len= *), intent( in) :: symbol_name

!  true if the symbol is a constant

logical, intent( in) :: is_const                                    

! **********************************************************************

!  entry: int_decl_str is blank_compress_lower_case integer declaration statement past the name
!         "" | "=..."
!         sym_name is the symbol name
!         is_const is true if this is a constant declaration
!         processing_set_file is true if this was found in the setfile

!  exit: integer declaration is added to the integer symbol list or error exit

! **********************************************************************

!  add_integer() local

! ----------------------------------------------------------------------

!  pointer to pre-existing symbols and the pointer to the new one

   type( symbol_t), pointer :: symbol_ptr

!  expression defining integer symbol

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  add_integer() text

continue

! ----------------------------------------------------------------------

!  check if already on integer list

   call seek_symbol_name( symbol_name, symbol_ptr)

   found_name: if( associated( symbol_ptr) )then
      
      name_type: select case( symbol_ptr% type_code)
   
      case( type_integer) name_type
   
         set_defn: if( symbol_ptr% predefined )then

            call msg_continue( "integer predeclared in setfile: " // trim( symbol_name) )

            param_match: if( is_const .neqv. symbol_ptr% constant )then

               call msg_quit( "inconsistent definition of symbol (constant v. variable): " &
                              // trim( symbol_name) )

            endif param_match

            symbol_ptr% predefined = .false.

            return

         else set_defn

            call msg_quit( "duplicate symbol declaration: " // trim( symbol_name) )

         endif set_defn

      case( type_logical) name_type
   
         call msg_quit( "integer name already defined as logical: " // trim( symbol_name) )

      case( type_macro) name_type
   
         call msg_quit( "integer name already defined as macro: " // trim( symbol_name) )

      case( type_text) name_type
   
         call msg_quit( "integer name already defined as text: " // trim( symbol_name) )

      end select name_type

   endif found_name

! ----------------------------------------------------------------------

!  insert a new symbol element into the symbol list

   call add_symbol( symbol_name, symbol_ptr)

! ----------------------------------------------------------------------

!  this symbol is an integer

   symbol_ptr% type_code = type_integer

!  store whether symbol is a parameter

   symbol_ptr% constant = is_const

!  store whether symbol is declared in the setfile

   symbol_ptr% predefined = processing_set_file

!  determine if declaration specifies a value

   got_eq: if( len( int_decl_str) > 0 )then

      symbol_ptr% defined = int_decl_str( : len( equals)) == equals

   else got_eq

      symbol_ptr% defined = .false.

   endif got_eq

! ----------------------------------------------------------------------

!  constants must have values

   constant_value: if( symbol_ptr% constant &
            .and. .not. symbol_ptr% defined )then

      call msg_quit( "an integer constant must have a value: " &
                     // trim( symbol_name) // trim( int_decl_str) )

   endif constant_value

!  decode the value

   process_value: if( symbol_ptr% defined )then

      all_constants = .true.

      expr_str = int_decl_str( len( equals) + 1: )

      call eval_int_expr( expr_str, symbol_ptr% integer_value)

      non_const: if( symbol_ptr% constant .and. .not. all_constants )then

         call msg_quit( "non constant expression used to define integer constant: " &
                        // trim( symbol_name))

      endif non_const

   endif process_value

! ----------------------------------------------------------------------

!  add_integer() exit

return

! **********************************************************************

!  add_integer()

end subroutine add_integer

! **********************************************************************
! **********************************************************************

!  add_symbol() allocate integer declaration in symbol table

subroutine add_symbol( symbol_name, symbol_ptr)

! **********************************************************************

!  add_symbol interface

! ----------------------------------------------------------------------

!  the name of the symbol to be added to the symbol list

character( len= *), intent( in) :: symbol_name

!  a pointer to the new entry

type( symbol_t), pointer :: symbol_ptr                              

! **********************************************************************

!  entry: none

!  exit: symbol_ptr points to a new element ready to be used

! **********************************************************************

!  add_symbol() local

! ----------------------------------------------------------------------

!  end of symbol list

   type( symbol_t), pointer, save :: current_symbol => null()

! ----------------------------------------------------------------------

!  catch allocate errors

   integer :: astat

! **********************************************************************

!  add_symbol() text

continue

! ----------------------------------------------------------------------

!  if already have symbol list

   start_or_append: if( associated( first_symbol) )then

      allocate( current_symbol% next, stat= astat)

      append_error: if( astat > 0 )then

         call msg_quit( "allocate symbol failed: " // trim( symbol_name))

      endif append_error

      current_symbol => current_symbol% next

   else start_or_append

      allocate( first_symbol, stat= astat)

      start_error: if( astat > 0 )then

         call msg_quit( "allocate symbol failed: " // trim( symbol_name))

      endif start_error

      current_symbol => first_symbol

   endif start_or_append

! ----------------------------------------------------------------------

!  update new entry's list components

   current_symbol% name_str = symbol_name

   current_symbol% type_code = type_none

   nullify( current_symbol% next)

!  initialize new entry's data components

   current_symbol% defined = .false.

   current_symbol% constant = .false.

   current_symbol% predefined = .false.

   current_symbol% integer_value = -huge( current_symbol% integer_value)

   current_symbol% logical_value = .false.

   nullify( current_symbol% dummy_args)

   current_symbol% macro_value = blank

   nullify( current_symbol% text_lines)

! ----------------------------------------------------------------------

!  set the caller's pointer to the new one

   symbol_ptr => current_symbol

! ----------------------------------------------------------------------

!  add_symbol() exit

return

! **********************************************************************

!  add_symbol()

end subroutine add_symbol

! **********************************************************************
! **********************************************************************

!  process_logical_declaration() process logical declarations

subroutine process_logical_declaration( logical_stmt)

! **********************************************************************

!  process_logical_declaration() interface

! ----------------------------------------------------------------------

!  the statement containing the logical declaration

character( len= *), intent( in) :: logical_stmt                    

! **********************************************************************

!  entry: logical_stmt is blank_compress_lower_case logical declaration past the logical keyword
!         "::..." | ",parameter::..."

!  exit: logical declaration is processed or error exit

! **********************************************************************

!  process_logical_declaration() local

! ----------------------------------------------------------------------

!  string containing a single symbol declaration symbol

   character( len= buffer_len) :: symbol_str

!  name of symbol

   character( len= symbol_name_len) :: symbol_name

!  results of decoding statement

   integer :: symbol_len

! ----------------------------------------------------------------------

!  point to next character to be decoded

   integer :: next_char

   integer :: decl_len

! **********************************************************************

!  process_logical_declaration() text

continue

   next_char = 1

! ----------------------------------------------------------------------

!  if active line, process the declaration

   active_line: if( if_construct% now_selected )then

!  extract all symbols on directive

      all_symbols: do

!  one symbol at a time to the symbol string

         decl_len = scan( logical_stmt( next_char: ), end_of_decl) + next_char - 1

         symbol_str = logical_stmt( next_char: decl_len - 1)

!  extract symbol name

         call get_symbol_name( symbol_str, symbol_name, symbol_len)

!  store symbol in symbol list

         call add_logical( symbol_str( symbol_len + 1: decl_len - 1), symbol_name, .false.)

!  comma separates symbols, blank is end of statement

         all_done: if( logical_stmt( decl_len: decl_len) == blank )then

            exit all_symbols

         endif all_done

!  reset for next symbol

         next_char = decl_len + 1

!  extract all symbols on directive

      enddo all_symbols

!  if active line, process the declaration

   endif active_line

! ----------------------------------------------------------------------

!  process_logical_declaration() exit

return

! **********************************************************************

!  process_logical_declaration()

end subroutine process_logical_declaration

! **********************************************************************
! **********************************************************************

!  process_logical_constant() process logical declarations

subroutine process_logical_constant( logical_stmt)

! **********************************************************************

!  process_logical_constant() interface

! ----------------------------------------------------------------------

!  the statement containing the declaration

character( len= *), intent( in) :: logical_stmt

! **********************************************************************

!  entry: logical_stmt is blank_compress_lower_case logical declaration past the logical keyword
!         "::..." | ",parameter::..."

!  exit: logical declaration is processed or error exit

! **********************************************************************

!  process_logical_constant() local

! ----------------------------------------------------------------------

!  string containing a single symbol declaration symbol

   character( len= buffer_len) :: symbol_str

!  name of symbol

   character( len= symbol_name_len) :: symbol_name

!  results of decoding statement

   integer :: symbol_len

! ----------------------------------------------------------------------

!  point to next character to be decoded

   integer :: next_char

   integer :: decl_len

! **********************************************************************

!  process_logical_constant() text

continue

   next_char = 1

! ----------------------------------------------------------------------

!  if active line, process the declaration

   active_line: if( if_construct% now_selected )then

!  extract all symbols on directive

      all_symbols: do

!  one symbol at a time to the symbol string

         decl_len = scan( logical_stmt( next_char: ), end_of_decl) + next_char - 1

         symbol_str = logical_stmt( next_char: decl_len - 1)

!  extract symbol name

         call get_symbol_name( symbol_str, symbol_name, symbol_len)

!  store symbol in symbol list

         call add_logical( symbol_str( symbol_len + 1: decl_len - 1), symbol_name, .true.)

!  comma separates symbols, blank is end of statement

         all_done: if( logical_stmt( decl_len: decl_len) == blank )then

            exit all_symbols

         endif all_done

!  reset for next symbol

         next_char = decl_len + 1

!  extract all symbols on directive

      enddo all_symbols

!  if active line, process the declaration

   endif active_line

! ----------------------------------------------------------------------

!  process_logical_constant() exit

return

! **********************************************************************

!  process_logical_constant()

end subroutine process_logical_constant

! **********************************************************************
! **********************************************************************

!  add_logical() copy integer symbol to symbol table

subroutine add_logical( log_decl_str, symbol_name, is_const)

! **********************************************************************

!  add_logical() interface

! ----------------------------------------------------------------------

!  the statement containing the logical declaration

character( len= *), intent( in) :: log_decl_str

!  the symbol name

character( len= *), intent( in) :: symbol_name

!  true if the symbol is a constant

logical, intent( in) :: is_const                                    

! **********************************************************************

!  entry: log_decl_str is blank_compress_lower_case logical declaration statement past the double colon
!         "name" | "name=..."
!         sym_name is the symbol name
!         is_const is true if this is a constant declaration
!         processing_set_file is true if this was found in the setfile

!  exit: logical declaration is added to the logical symbol list or error exit

! **********************************************************************

!  add_logical() local

! ----------------------------------------------------------------------

!  pointer to pre-existing symbol

   type( symbol_t), pointer :: symbol_ptr

!  evaluate expression if there is one

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  add_logical() text

continue

! ----------------------------------------------------------------------

!  check if already on logical list

   call seek_symbol_name( symbol_name, symbol_ptr)

   duplicate_log: if( associated( symbol_ptr) )then

      name_type: select case( symbol_ptr% type_code)
   
      case( type_logical) name_type
   
         set_defn: if( symbol_ptr% predefined )then

            call msg_continue( "logical predeclared in setfile: " // trim( symbol_name) )

            param_match: if( is_const .neqv. symbol_ptr% constant )then

               call msg_quit( "inconsistent definition of logical (constant v. variable): " &
                              // trim( symbol_name) )

            endif param_match

            symbol_ptr% predefined = .false.

            return

         else set_defn

            call msg_quit( "duplicate logical declaration: " // trim( symbol_name) )

         endif set_defn

      case( type_integer) name_type
   
         call msg_quit( "logical name already defined as integer: " // trim( symbol_name) )

      case( type_macro) name_type
   
         call msg_quit( "integer name already defined as macro: " // trim( symbol_name) )

      case( type_text) name_type
   
         call msg_quit( "integer name already defined as text: " // trim( symbol_name) )

      end select name_type

   endif duplicate_log

! ----------------------------------------------------------------------

!  allocate next logical

   call add_symbol( symbol_name, symbol_ptr)

! ----------------------------------------------------------------------

!  store symbol code

   symbol_ptr% type_code = type_logical

!  store whether symbol is a parameter

   symbol_ptr% constant = is_const

!  store whether symbol is declared in the setfile

   symbol_ptr% predefined = processing_set_file

!  determine if declaration specifies a value

   got_eq: if( len( log_decl_str) > 0 )then

      symbol_ptr% defined = log_decl_str( : len( equals)) == equals

   else got_eq

      symbol_ptr% defined = .false.

   endif got_eq

! ----------------------------------------------------------------------

!  constants must have values

   constant_value: if( symbol_ptr% constant &
            .and. .not. symbol_ptr% defined )then

         call msg_quit( "a logical constant must have a value: " &
                        // trim( symbol_name) // trim( log_decl_str) )

   endif constant_value

!  decode the value

   got_value: if( symbol_ptr% defined )then

      all_constants = .true.

      expr_str = log_decl_str( len( equals) + 1: )

      call eval_log_expr( expr_str, symbol_ptr% logical_value)

!  evaluate the expression

      non_const: if( symbol_ptr% constant .and. .not. all_constants )then

         call msg_quit( "non constant expression used to define logical constant: " &
                        // trim( symbol_name) // trim( expr_str))

      endif non_const

   endif got_value

! ----------------------------------------------------------------------

!  add_logical() exit

return

! **********************************************************************

!  add_logical()

end subroutine add_logical

! **********************************************************************
! **********************************************************************

!  %%% process coco constructs

! **********************************************************************
! **********************************************************************

!  process_macro_declaration() process macro declarations

subroutine process_macro_declaration( mac_stmt)

! **********************************************************************

!  process_macro_declaration() interface

! ----------------------------------------------------------------------

!  the statement containing the macro declaration

character( len= *), intent( in) :: mac_stmt                        

! **********************************************************************

!  entry: mac_stmt is blank_compress_lower_case logical declaration
!         past the macro keyword "::..."

!  exit: macro declaration is processed or error exit

! **********************************************************************

!  process_macro_declaration() local

! ----------------------------------------------------------------------

!  name of symbol

   character( len= symbol_name_len) :: symbol_name

!  results of decoding statement

   integer :: symbol_len

! ----------------------------------------------------------------------

!  point to next character to be decoded

   integer :: next_char

! **********************************************************************

!  process_macro_declaration() text

continue

   next_char = 1

! ----------------------------------------------------------------------

!  if active line, process the declaration

   active_line: if( if_construct% now_selected )then

!  extract symbol name

      call get_macro_name( mac_stmt( next_char: ), symbol_name, symbol_len)

!  store symbol in symbol list

      next_char = next_char + symbol_len

      call add_macro( mac_stmt( next_char: ), symbol_name)

! ----------------------------------------------------------------------

!  if reporting extensions

      extensions: if( options% report_extensions )then

         call msg_continue( "processed macro directive: " // trim( mac_stmt))

      endif extensions

   endif active_line

! ----------------------------------------------------------------------

!  process_macro_declaration() exit

return

! **********************************************************************

!  process_macro_declaration()

end subroutine process_macro_declaration

! **********************************************************************
! **********************************************************************

!  process_undefine_directive() process macro declarations

subroutine process_undefine_directive( undefine_stmt)

! **********************************************************************

!  process_undefine_directive() interface

! ----------------------------------------------------------------------

!  the statement containing the macro declaration

character( len= *), intent( in) :: undefine_stmt                        

! **********************************************************************

!  entry: undefine_stmt is blank_compress_lower_case logical declaration
!         past the undefine: keyword: "symbol[,symbol]..."

!  exit: macro declaration is processed or error exit

! **********************************************************************

!  process_undefine_directive() local

! ----------------------------------------------------------------------

!  name of symbol

   character( len= symbol_name_len) :: symbol_name

!  results of decoding statement

   integer :: symbol_len

! ----------------------------------------------------------------------

!  pointers to characters

   integer :: next_char

! **********************************************************************

!  process_undefine_directive() text

continue

! ----------------------------------------------------------------------

!  syntax check

   not_null: if( undefine_stmt == blank )then

      call msg_quit( "no names on undefine statement")

   endif not_null

! ----------------------------------------------------------------------

!  if active line, process the declaration

   active_line: if( if_construct% now_selected )then

!  extract all symbols on directive

      next_char = 1

      all_symbols: do

!  one symbol at a time to the symbol string

         symbol_len = scan( undefine_stmt( next_char: ), end_of_decl) + next_char - 1

         symbol_name = undefine_stmt( next_char: symbol_len - 1)

!  store symbol in symbol list

         call remove_symbol( symbol_name)

!  comma separates symbols, blank is end of statement

         all_done: if( undefine_stmt( symbol_len: symbol_len) == blank )then

            exit all_symbols

         endif all_done

!  reset for next symbol

         next_char = symbol_len + 1

!  extract all symbols on directive

      enddo all_symbols

! ----------------------------------------------------------------------

!  if reporting extensions

      extensions: if( options% report_extensions )then

         call msg_continue( "processed undefine directive: " // trim( undefine_stmt))

      endif extensions

   endif active_line

! ----------------------------------------------------------------------

!  process_undefine_directive() exit

return

! **********************************************************************

!  process_undefine_directive()

end subroutine process_undefine_directive

! **********************************************************************
! **********************************************************************

!  remove_symbol() seek symbol on symbol list

subroutine remove_symbol( symbol_str)

! **********************************************************************

!  remove_symbol() interface

! ----------------------------------------------------------------------

!  the name of the symbol being sought

character( len= *), intent( in) :: symbol_str

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  remove_symbol() local

! ----------------------------------------------------------------------

!  a pointer to the symbol found

   type( symbol_t), pointer :: symbol_ptr

!  preserve the previous pointer

   type( symbol_t), pointer :: previous_ptr

! ----------------------------------------------------------------------

!  check deallocation status

   integer :: astat

! **********************************************************************

!  remove_symbol() text

continue

! ----------------------------------------------------------------------

!  search symbol list

   nullify( symbol_ptr)

   nullify( previous_ptr)

   symbol_ptr => first_symbol

!  loop through symbol list

   search_list: do while( associated( symbol_ptr) )

!  if found name matching

      name_match: if( symbol_str == symbol_ptr% name_str )then

         relink_list: if( associated( previous_ptr) )then

            previous_ptr% next => symbol_ptr% next

         else relink_list

            first_symbol => symbol_ptr% next

         endif relink_list

!  deallocate any allocated components

         dealloc_args: if( associated( symbol_ptr% dummy_args) )then

            deallocate( symbol_ptr% dummy_args, stat= astat)

            dummy_args_error: if( astat > 0 )then

               call msg_quit( "error deallocating dummy arguments: " // trim( symbol_ptr% dummy_args( 1)))

            endif dummy_args_error

         endif dealloc_args

         dealloc_text: if( associated( symbol_ptr% text_lines) )then

            deallocate( symbol_ptr% text_lines, stat= astat)

            text_lines_error: if( astat > 0 )then

               call msg_quit( "error deallocating text block: " // trim( symbol_ptr% text_lines( 1)))

            endif text_lines_error

         endif dealloc_text

!  deallocate its entry

         deallocate( symbol_ptr, stat= astat)

         dealloc_error: if( astat > 0 )then

            call msg_quit( "error deallocating symbol: " // trim( symbol_str))

         endif dealloc_error

         return

      endif name_match

!  else update the pointers and try again

      update_ptr: if( associated( previous_ptr) )then

         previous_ptr => symbol_ptr

         symbol_ptr => symbol_ptr% next

      else update_ptr

         previous_ptr => first_symbol

         symbol_ptr => first_symbol% next

      endif update_ptr

   enddo search_list

!  searched entire list with no match for name   

   call msg_quit( "symbol not found: " // trim( symbol_str))

! ----------------------------------------------------------------------

!  remove_symbol() exit

return

! **********************************************************************

!  remove_symbol()

end subroutine remove_symbol

! **********************************************************************
! **********************************************************************

!  add_macro() copy integer symbol to symbol table

subroutine add_macro( mac_decl_str, symbol_name)

! **********************************************************************

!  add_macro() interface

! ----------------------------------------------------------------------

!  the statement containing the macro definition

character( len= *), intent( in) :: mac_decl_str

!  the macro name

character( len= *), intent( in) :: symbol_name                    

! **********************************************************************

!  entry: mac_decl_str is blank_compress_lower_case logical declaration statement past the double colon
!         "" | "=..."
!         sym_name is the symbol name
!         is_const is true if this is a constant declaration
!         processing_set_file is true if this was found in the setfile

!  exit: logical declaration is added to the logical symbol list or error exit

! **********************************************************************

!  add_macro() local

! ----------------------------------------------------------------------

!  pointer to pre-existing symbol and to the new symbol

   type( symbol_t), pointer :: symbol_ptr

!  location of close parenthesis

   integer :: close_idx

   integer :: char_idx

!  buffer for macro argument list

   character( len= buffer_len) :: arg_list_str

! **********************************************************************

!  add_macro() text

continue

! ----------------------------------------------------------------------

!  check if already on symbol list

   call seek_symbol_name( symbol_name, symbol_ptr)

   duplicate_mac: if( associated( symbol_ptr) )then

      name_type: select case( symbol_ptr% type_code)
   
      case( type_macro) name_type
   
         call msg_quit( "duplicate macro declaration: " // trim( symbol_name) )

      case( type_integer) name_type
   
         call msg_quit( "macro name already defined as integer: " // trim( symbol_name) )

      case( type_logical) name_type
   
         call msg_quit( "macro name already defined as macro: " // trim( symbol_name) )

      case( type_text) name_type
   
         call msg_quit( "macro name already defined as text: " // trim( symbol_name) )

      end select name_type

   endif duplicate_mac

! ----------------------------------------------------------------------

!  allocate next macro

   call add_symbol( symbol_name, symbol_ptr)

! ----------------------------------------------------------------------

!  store symbol code

   symbol_ptr% type_code = type_macro

!  determine if declaration specifies an arglist

   got_arglist: if( mac_decl_str( len( open_paren): len( open_paren)) == open_paren )then

! paren enclose arglist

      call seek_close_paren( mac_decl_str, len( open_paren), close_idx)

      no_close_paren: if( close_idx == 0 )then

         call msg_quit( "no closing parenthesis: " // trim( mac_decl_str) )

      endif no_close_paren

      arg_list_str = mac_decl_str( len( open_paren) + 1: close_idx - 1)

      call process_dummy_arglist( arg_list_str, symbol_ptr)

      call verify_dummy_args( symbol_ptr)

      char_idx = close_idx + len( equals) + 1

   else got_arglist

      char_idx = len( equals) + 1

      nullify( symbol_ptr% dummy_args)

   endif got_arglist

!  validate macro definition

   call verify_macro_value( symbol_ptr, mac_decl_str( char_idx: ))

!  store the value

   symbol_ptr% macro_value = mac_decl_str( char_idx: )

! ----------------------------------------------------------------------

!  add_macro() exit

return

! **********************************************************************

!  add_macro()

end subroutine add_macro

! **********************************************************************
! **********************************************************************

!  process_dummy_arglist() process macro or text dummy argument list

subroutine process_dummy_arglist( arglist, symbol_ptr)

! **********************************************************************

!  process_dummy_arglist() interface

! ----------------------------------------------------------------------

!  the comma separated dummy argument list

character( len= *), intent( in) :: arglist

!  a pointer to the macro

type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  entry: arglist is a character with the arglist

!  exit: symbol_ptr has its arglist array allocated or defined

! **********************************************************************

!  process_dummy_arglist() constants

   character( len= *), parameter :: end_of_arg = blank // comma

! **********************************************************************

!  process_dummy_arglist() local

! ----------------------------------------------------------------------

   character( len= buffer_len) :: arg_str

!  length of dummy argument list

   integer :: arg_len

! ----------------------------------------------------------------------

!  number of dummy arguments found

   integer :: number_of_args

!  point to next character to be decoded

   integer :: next_char

   integer :: this_arg

!  allocation status

   integer :: astat

! **********************************************************************

!  process_dummy_arglist() text

continue

! ----------------------------------------------------------------------

!  get length of arglist

   arg_str = arglist

   arg_len = len_trim( arglist)

!  if the arglist has characters

   null_list: if( arg_len > 0 )then

!  count macro dummy arguments

      number_of_args = 1

      count_args: do next_char = 1, arg_len

         found_another: select case( arg_str( next_char: next_char))

         case( comma) found_another

            number_of_args = number_of_args + 1

         end select found_another

      enddo count_args

!  allocate array to hold dummy args

      allocate( symbol_ptr% dummy_args( number_of_args), stat= astat)

      arg_error: if( astat > 0 )then

         call msg_quit( "allocate dummy argument list failed: " // trim( symbol_ptr% name_str))

      endif arg_error

!  put each dummy arg into its own array element

      get_args: do this_arg = 1, number_of_args

!  find end of each arg

         next_char = scan( arg_str, end_of_arg)

!  store arg

         symbol_ptr% dummy_args( this_arg) = arg_str( : next_char - 1)

!  shift off that arg

         arg_str = arg_str( next_char + 1: )

      enddo get_args

   else null_list

      nullify( symbol_ptr% dummy_args)

   endif null_list

! ----------------------------------------------------------------------

!  process_dummy_arglist() exit

return

! **********************************************************************

!  process_dummy_arglist()

end subroutine process_dummy_arglist

! **********************************************************************
! **********************************************************************

!  verify_macro_value() process macro dummy argument list

subroutine verify_macro_value( macro_ptr, value_str)

! **********************************************************************

!  verify_macro_value() interface

! ----------------------------------------------------------------------

!  a pointer to the macro

type( symbol_t), pointer :: macro_ptr

!  the macro value string

character( len=*), intent( in) :: value_str                       

! **********************************************************************

!  entry: macro or its coridate value

!  exit: macro's value is valid & may be assigned to macro

! **********************************************************************

!  verify_macro_value() local

! ----------------------------------------------------------------------

!  make ?name? string

   character( len= target_len) :: search_str

!  point to characters

   integer :: dup_idx

!  point to dummy args

   integer :: this_arg

! **********************************************************************

!  verify_macro_value() text

continue

! ----------------------------------------------------------------------

!  check for null valued macros

   null_macro: if( len_trim( value_str) == 0 )then

      call msg_quit( "macro has null value: " // trim( macro_ptr% name_str) )

   endif null_macro

! ----------------------------------------------------------------------

!  check for recursive macros

   search_str = arg_key // trim( macro_ptr% name_str) // arg_key

   dup_idx = index( value_str, trim( search_str))

   recursive_mac: if( dup_idx > 0 )then

      call msg_quit( "recursive macro definition: " // trim( macro_ptr% name_str) )

   endif recursive_mac

!  check that dummy args all appear in the macro value

   have_dummy_args: if( associated( macro_ptr% dummy_args) )then

      scan_dummy_args: do this_arg = lbound( macro_ptr% dummy_args, dim= 1), ubound( macro_ptr% dummy_args, dim= 1)

         search_str = arg_key // trim( macro_ptr% dummy_args( this_arg)) // arg_key

         dup_idx = index( value_str, trim( search_str))

         arg_not_used: if( dup_idx == 0 )then

            call msg_quit( "macro argument unused in: " // trim( macro_ptr% name_str) &
                           // " argument: " // trim( macro_ptr% dummy_args( this_arg)) )

         endif arg_not_used

      enddo scan_dummy_args

   endif have_dummy_args

! ----------------------------------------------------------------------

!  verify_macro_value() exit

return

! **********************************************************************

!  verify_macro_value()

end subroutine verify_macro_value

! **********************************************************************
! **********************************************************************

!  verify_dummy_args() process macro or text dummy argument list

subroutine verify_dummy_args( symbol_ptr)

! **********************************************************************

!  verify_dummy_args() interface

! ----------------------------------------------------------------------

!  a pointer to the macro

type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  entry: macro or text dummy arg list

!  exit: macro's value is valid & may be assigned to macro

! **********************************************************************

!  verify_dummy_args() local

! ----------------------------------------------------------------------

!  point to characters

   integer :: this_arg

   integer :: dup_idx

   integer :: char_idx

   integer :: arg_name_len

! ----------------------------------------------------------------------

!  check for the same name amongst integers, logicals, macros, text blocks

   type( symbol_t), pointer :: name_ptr

! **********************************************************************

!  verify_dummy_args() text

continue

! ----------------------------------------------------------------------

!  macro or text dummy arguments must be valid names

   check_names: do this_arg = lbound( symbol_ptr% dummy_args, dim= 1), ubound( symbol_ptr% dummy_args, dim= 1)

!  check whether initial character is alphabetic

      char_idx = verify( symbol_ptr% dummy_args( this_arg)( 1: 1), alpha_chars)

      initial_ok: if( char_idx /= 0 )then

         call msg_quit( "illegal initial character in macro arg name: " // &
                        trim( symbol_ptr% dummy_args( this_arg)) )

      endif initial_ok

!  check whether following characters in name are legal

      arg_name_len = len_trim( symbol_ptr% dummy_args( this_arg))

      char_idx = verify( symbol_ptr% dummy_args( this_arg)( 2: arg_name_len), alphanum_chars)

      name_ok: if( char_idx /= 0 )then

         call msg_quit( "illegal character in macro arg name: " // &
                        trim( symbol_ptr% dummy_args( this_arg)) )

      endif name_ok

!  check whether names of args duplicates arg names

      self_dup: do dup_idx = 1, this_arg - 1

         dup_arg: if( symbol_ptr% dummy_args( this_arg) == symbol_ptr% dummy_args( dup_idx) )then

            call msg_quit( "duplicate dummy argument name: " // &
                           symbol_ptr% dummy_args( this_arg) )

         endif dup_arg

      enddo self_dup

!  check whether arg name duplicates symbol names

      call seek_symbol_name( symbol_ptr% dummy_args( this_arg), name_ptr)

      dup_sym: if( associated( name_ptr) )then

         call msg_quit( "dummy argument name duplicates symbol name: " // &
                        symbol_ptr% dummy_args( this_arg) )

      endif dup_sym

   enddo check_names

! ----------------------------------------------------------------------

!  verify_dummy_args() exit

return

! **********************************************************************

!  verify_dummy_args()

end subroutine verify_dummy_args

! **********************************************************************
! **********************************************************************

!  %%% process if constructs: if, elseif, else, endif, ifndef

! **********************************************************************
! **********************************************************************

!  process_if_directive() process a coco if( )then directive

subroutine process_if_directive( if_dir)

! **********************************************************************

!  process_if_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an if directive

character( len= *), intent( in) :: if_dir                          

! **********************************************************************

!  entry: if_dir is blank_compress_lower_case if directive, past the coco key and the "if("
!         "<logical>)then"

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_if_directive() local

! ----------------------------------------------------------------------

!  pointer to ')then'

   integer :: then_idx

!  value of logical expression

   logical :: expression_value

!  status of allocating a new if block

   integer :: astat

!  copy expression string for evaluation

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  process_if_directive() text

continue

! ----------------------------------------------------------------------

!  check for a well formed if directive

   then_idx = index( if_dir, then_str)

   syntax_check: if( then_idx == 0 )then

      call msg_quit( "no 'then' in if directive: " // trim( if_dir) )

   endif syntax_check

   extra_chars_check: if( if_dir( then_idx + len( then_str): ) /= blank )then

      call msg_quit( "extra characters after if directive: " // trim( if_dir) )

   endif extra_chars_check

! ----------------------------------------------------------------------

!  append new if construct at end of list

   allocate( if_construct% nested, stat= astat)

   alloc_error: if( astat > 0 )then

      call msg_quit( "allocate if block failed: " // trim( if_dir) )

   endif alloc_error

!  establish pointers

   if_construct% nested% enclosing => if_construct

!  make new if construct the active if construct

   if_construct => if_construct% nested

   nullify( if_construct% nested)

!  this phase is an if block

   if_construct% phase = if_block

! ----------------------------------------------------------------------

!  if this if block is enclosed within selected lines

   active_lines: if( if_construct% enclosing% now_selected )then

! ----------------------------------------------------------------------

!  evaluate logical expression only when enclosing if block is selected

      expr_str = if_dir( : then_idx - 1)

      call eval_log_expr( expr_str, expression_value)

!  set if value accordingly

      if_construct% now_selected = expression_value

      if_construct% ever_selected = expression_value

! ----------------------------------------------------------------------

!  the enclosing if block is not selected

   else active_lines

      if_construct% now_selected = .false.

      if_construct% ever_selected = .true.

   endif active_lines

! ----------------------------------------------------------------------

!  process_if_directive() exit

return

! **********************************************************************

!  process_if_directive()

end subroutine process_if_directive

! **********************************************************************
! **********************************************************************

!  process_elseif_directive() process a coco elseif( )then directive

subroutine process_elseif_directive( elseif_dir)

! **********************************************************************

!  process_elseif_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an elseif directive

character( len= *), intent( in) :: elseif_dir                      

! **********************************************************************

!  entry: elseif_dir is blank_compress_lower_case elseif directive, past the coco key

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_elseif_directive() local

! ----------------------------------------------------------------------

!  location of closing )then

   integer :: then_idx

!  value of logical expression

   logical :: expression_value

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  process_elseif_directive() text

continue

! ----------------------------------------------------------------------

!  if not in if-block, elseif is misplaced

   if_sequence: select case( if_construct% phase)

   case( outside_block ) if_sequence

      call msg_quit( "elseif outside if construct: " // trim( elseif_dir) )

   case( else_block ) if_sequence

      call msg_quit( "elseif after else: " // trim( elseif_dir) )

   end select if_sequence

! ----------------------------------------------------------------------

!  logical expression must be between 'if(' and ')then'

   then_idx = index( elseif_dir, then_str)

   syntax_check: if( then_idx == 0 )then

      call msg_quit( "no 'then' in elseif directive: " // trim( elseif_dir) )

   endif syntax_check

   extra_chars_check: if( elseif_dir( then_idx + len( then_str): ) /= blank )then

      call msg_quit( "extra characters after elseif directive: " // trim( elseif_dir) )

   endif extra_chars_check

!  this phase is an elseif block

   if_construct% phase = elseif_block

! ----------------------------------------------------------------------

!  if this if block is enclosed within selected lines

   active_lines: if( if_construct% enclosing% now_selected )then

      previous_true: if( if_construct% ever_selected )then

         if_construct% now_selected = .false.

      else previous_true

!  evaluate logical expression

         expr_str = elseif_dir( : then_idx - 1)

         call eval_log_expr( expr_str, expression_value)

!  set if value accordingly

         if_construct% now_selected = expression_value

         if_construct% ever_selected = expression_value

      endif previous_true

   endif active_lines

! ----------------------------------------------------------------------

!  process_elseif_directive() exit

return

! **********************************************************************

!  process_elseif_directive()

end subroutine process_elseif_directive

! **********************************************************************
! **********************************************************************

!  process_else_directive() process a coco else directive

subroutine process_else_directive( else_dir)

! **********************************************************************

!  process_else_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an else directive

character( len= *), intent( in) :: else_dir                         

! **********************************************************************

!  entry: else_dir is blank_compress_lower_case else directive, past the coco key

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_else_directive() text

continue

! ----------------------------------------------------------------------

!  if not in if-block, else is misplaced

   if_sequence: select case( if_construct% phase)

   case( outside_block ) if_sequence

      call msg_quit( "else outside if construct: " // trim( else_dir) )

   case( else_block ) if_sequence

      call msg_quit( "else after else: " // trim( else_dir) )

   end select if_sequence

!  must have nothing after 'else'

   syntax_error: if( else_dir /= blank )then

      call msg_quit( "extra characters after else directive: " // trim( else_dir) )

   endif syntax_error

! ----------------------------------------------------------------------

!  this phase is an else block

   if_construct% phase = else_block

! ----------------------------------------------------------------------

!  select else block if this if ithe enclosing block is active and no previous block has been selected

   if_construct% now_selected = if_construct% enclosing% now_selected &
                                .and. .not. if_construct% ever_selected

! ----------------------------------------------------------------------

!  process_else_directive() exit

return

! **********************************************************************

!  process_else_directive()

end subroutine process_else_directive

! **********************************************************************
! **********************************************************************

!  process_endif_directive() process a coco endif directive

subroutine process_endif_directive( endif_dir)

! **********************************************************************

!  process_endif_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an endif directive

character( len= *), intent( in) :: endif_dir                     

! **********************************************************************

!  entry: endif_dir is blank_compress_lower_case endif directive, past the coco key

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_endif_directive() local

! ----------------------------------------------------------------------

!  deallocate status

   integer :: astat

! **********************************************************************

!  process_endif_directive() text

continue

! ----------------------------------------------------------------------

!  if not in if-block, endif is misplaced

   if_sequence: select case( if_construct% phase)

   case( outside_block ) if_sequence

      call msg_quit( "endif outside any if construct: " // trim( endif_dir) )

   end select if_sequence

!  must have nothing after 'endif'

   syntax_error: if( endif_dir /= blank )then

      call msg_quit( "extra characters after endif directive: " // trim( endif_dir) )

   endif syntax_error

! ----------------------------------------------------------------------

!  decrement if level

   if_construct => if_construct% enclosing

   deallocate( if_construct% nested, stat= astat)

   next_error: if( astat > 0 )then

      call msg_quit( "deallocate if construct failed: " // trim( endif_dir) )

   endif next_error

   nullify( if_construct% nested)

! ----------------------------------------------------------------------

!  process_endif_directive() exit

return

! **********************************************************************

!  process_endif_directive()

end subroutine process_endif_directive

! **********************************************************************
! **********************************************************************

!  process_ifdef_directive() process a coco ifdef( )then directive

subroutine process_ifdef_directive( ifdef_dir)

! **********************************************************************

!  process_ifdef_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an if directive

character( len= *), intent( in) :: ifdef_dir

! **********************************************************************

!  entry: ifdef_dir is blank_compress_lower_case ifdef directive,
!         past the coco key and the "ifdef("
!         "<logical>)then"

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_ifdef_directive() local

! ----------------------------------------------------------------------

!  pointer to ')then'

   integer :: then_idx

!  status of allocating a new if block

   integer :: astat

!  copy symbol name string to see if it exists

   character( len= symbol_name_len) :: name_str

!  symbol pointer to return symbol or null

   type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  process_ifdef_directive() text

continue

! ----------------------------------------------------------------------

!  check for a well formed if directive

   then_idx = index( ifdef_dir, then_str)

   syntax_check: if( then_idx == 0 )then

      call msg_quit( "no 'then' in ifdef directive: " // trim( ifdef_dir) )

   endif syntax_check

   extra_chars_check: if( ifdef_dir( then_idx + len( then_str): ) /= blank )then

      call msg_quit( "extra characters after ifdef directive: " // trim( ifdef_dir) )

   endif extra_chars_check

! ----------------------------------------------------------------------

!  append new if construct at end of list

   allocate( if_construct% nested, stat= astat)

   alloc_error: if( astat > 0 )then

      call msg_quit( "allocate ifdef block failed: " // trim( ifdef_dir) )

   endif alloc_error

!  establish pointers

   if_construct% nested% enclosing => if_construct

!  make new if construct the active if construct

   if_construct => if_construct% nested

   nullify( if_construct% nested)

!  this phase is an if block

   if_construct% phase = else_block

! ----------------------------------------------------------------------

!  if this if block is enclosed within selected lines

   active_lines: if( if_construct% enclosing% now_selected )then

! ----------------------------------------------------------------------

!  check existance only when enclosing if block is selected

      name_str = ifdef_dir( : then_idx - 1)

      call seek_symbol_name( name_str, symbol_ptr)

!  set ifdef value accordingly

      if_construct% now_selected = associated( symbol_ptr)

      if_construct% ever_selected = associated( symbol_ptr)

! ----------------------------------------------------------------------

!  the enclosing if block is not selected

   else active_lines

      if_construct% now_selected = .false.

      if_construct% ever_selected = .true.

   endif active_lines

! ----------------------------------------------------------------------

!  process_ifdef_directive() exit

return

! **********************************************************************

!  process_ifdef_directive()

end subroutine process_ifdef_directive

! **********************************************************************
! **********************************************************************

!  process_ifndef_directive() process a coco ifndef( )then directive

subroutine process_ifndef_directive( ifndef_dir)

! **********************************************************************

!  process_ifndef_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an if directive

character( len= *), intent( in) :: ifndef_dir

! **********************************************************************

!  entry: ifndef_dir is blank_compress_lower_case ifndef directive,
!         past the coco key and the "ifndef("
!         "<logical>)then"

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_ifndef_directive() local

! ----------------------------------------------------------------------

!  pointer to ')then'

   integer :: then_idx

!  status of allocating a new if block

   integer :: astat

!  copy symbol name string to see if it exists

   character( len= symbol_name_len) :: name_str

!  symbol pointer to return symbol or null

   type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  process_ifndef_directive() text

continue

! ----------------------------------------------------------------------

!  check for a well formed if directive

   then_idx = index( ifndef_dir, then_str)

   syntax_check: if( then_idx == 0 )then

      call msg_quit( "no 'then' in ifndef directive: " // trim( ifndef_dir) )

   endif syntax_check

   extra_chars_check: if( ifndef_dir( then_idx + len( then_str): ) /= blank )then

      call msg_quit( "extra characters after ifndef directive: " // trim( ifndef_dir) )

   endif extra_chars_check

! ----------------------------------------------------------------------

!  append new if construct at end of list

   allocate( if_construct% nested, stat= astat)

   alloc_error: if( astat > 0 )then

      call msg_quit( "allocate ifndef block failed: " // trim( ifndef_dir) )

   endif alloc_error

!  establish pointers

   if_construct% nested% enclosing => if_construct

!  make new if construct the active if construct

   if_construct => if_construct% nested

   nullify( if_construct% nested)

!  this phase is an if block

   if_construct% phase = else_block

! ----------------------------------------------------------------------

!  if this if block is enclosed within selected lines

   active_lines: if( if_construct% enclosing% now_selected )then

! ----------------------------------------------------------------------

!  check existance only when enclosing if block is selected

      name_str = ifndef_dir( : then_idx - 1)

      call seek_symbol_name( name_str, symbol_ptr)

!  set ifndef value accordingly

      if_construct% now_selected = .not. associated( symbol_ptr)

      if_construct% ever_selected = .not. associated( symbol_ptr)

! ----------------------------------------------------------------------

!  the enclosing if block is not selected

   else active_lines

      if_construct% now_selected = .false.

      if_construct% ever_selected = .true.

   endif active_lines

! ----------------------------------------------------------------------

!  process_ifndef_directive() exit

return

! **********************************************************************

!  process_ifndef_directive()

end subroutine process_ifndef_directive

! **********************************************************************
! **********************************************************************

!  process_assert_directive() process an assert directive

subroutine process_assert_directive( assert_dir)

! **********************************************************************

!  process_assert_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an assert directive

character( len= *), intent( in) :: assert_dir

! **********************************************************************

!  entry: assert_dir is quoted assert condition

!  exit: assertion code is written to the output

! **********************************************************************

!  process_assert_directive() constants

! ----------------------------------------------------------------------

!  pieces of the assert output

   character( len= *), parameter :: if_prt = 'if( .not. ( '

   character( len= *), parameter :: then_prt = ' ) )then'

   character( len= *), parameter :: write_prt = 'write( unit= error_unit, fmt= *) "assertion failed: '

   character( len= *), parameter :: stop_prt = 'stop "assertion failed"'

   character( len= *), parameter :: endif_prt = 'endif'

! ----------------------------------------------------------------------

!  starting column of output

   integer, parameter :: start_col = 7

! **********************************************************************

!  process_assert_directive() local

! ----------------------------------------------------------------------

!  the condition must be unquoted before use

   character( len= buffer_len) :: assert_cond

!  assemble the output line

   character( len= buffer_len) :: edit_line

! ----------------------------------------------------------------------

!  convert the current input line number to characters

   character( len= conversion_len) :: conversion_str

!  lengths of quoted and unquoted include file name

   integer :: construct_len

   integer :: cond_len

! ----------------------------------------------------------------------

!  inhibit wrapping when editing, assert processing doesn its own wrapping

   integer :: wrap_state

! **********************************************************************

!  process_assert_directive() text

continue

! ----------------------------------------------------------------------

!  check syntax- condition must be string within quotes

   call unquote_string( assert_dir, assert_cond, construct_len, cond_len )

   no_cond_str: if( cond_len == 0 .or. construct_len == 0 )then

      call msg_quit( "can't find assert condition: " // trim( assert_dir))

   endif no_cond_str

!  check syntax- directive must be blank after condition

   extra_chars: if( assert_dir( construct_len + 1: ) /= blank )then

      call msg_quit( "extra characters after assert condition: " // trim( assert_dir))

   endif extra_chars

! ----------------------------------------------------------------------

!  if active block, process assert directive

   active_line: if( if_construct% now_selected )then

!  editing source lines is enabled

      edit_source: if( options% edit_source )then

!  if ? is present edit source line

         edit_cond: if( index( assert_cond, arg_key) > 0 )then

!  disable wrapping because each line will be wrapped after it is constructed

            wrap_state = options% wrap_length

            options% wrap_length = wrap_off

            call edit_source_line( assert_cond)

            options% wrap_length = wrap_state

         endif edit_cond

      endif edit_source

! ----------------------------------------------------------------------

!  write the if statement

      edit_line( : start_col - 1) = blank

      edit_line( start_col: ) = if_prt // trim( assert_cond) // then_prt

!  remove any line length overflow

      wrap_if: if( options% wrap_length /= wrap_off )then

         call wrap_source_line( edit_line)

      endif wrap_if

! ----------------------------------------------------------------------

!  write assembled if-then statement

      line = edit_line

      call write_source_line( output_file)

! ----------------------------------------------------------------------

!  write the write statement

      edit_line( : start_col - 1) = blank

!  get the current line number

      write( unit= conversion_str, fmt= conversion_fmt) current_file% lines_transfered

!  construct the assertion complaint

      edit_line( start_col: ) = write_prt // trim( current_file% name_str) &
                                          // ": " // trim( adjustl( conversion_str)) &
                                          // ': " // ' // trim( assert_dir)

!  remove any line length overflow

      wrap_write: if( options% wrap_length /= wrap_off )then

         call wrap_source_line( edit_line)

      endif wrap_write

! ----------------------------------------------------------------------

!  write assembled write statement

      line = edit_line

      call write_source_line( output_file)

! ----------------------------------------------------------------------

!  blank until the start column for the stop and endif

      line( : start_col - 1) = blank

!  write the stop statement

      line( start_col: ) = stop_prt

      call write_source_line( output_file)

! ----------------------------------------------------------------------

!  write the endif statement

      line( start_col: ) = endif_prt

      call write_source_line( output_file)

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed assert directive")

      endif extensions

! ----------------------------------------------------------------------

   endif active_line

! ----------------------------------------------------------------------

!  process_assert_directive() exit

return

! **********************************************************************

!  process_assert_directive()

end subroutine process_assert_directive

! **********************************************************************
! **********************************************************************

!  process_assertif_directive() process an assertif directive

subroutine process_assertif_directive( assertif_dir)

! **********************************************************************

!  process_assertif_directive() interface

! ----------------------------------------------------------------------

!  a statement containing an assertif directive

character( len= *), intent( in) :: assertif_dir

! **********************************************************************

!  entry: assertif_dir "<logical expression>)'assert condition'"

!  exit: if logical expression, assertion code is written to the output

! **********************************************************************

!  process_assertif_directive() local

! ----------------------------------------------------------------------

!  length of logical condition

   integer :: logical_len

!  value of logical condition

   logical :: logical_value

! **********************************************************************

!  process_assertif_directive() text

continue

! ----------------------------------------------------------------------

!  check syntax- condition must be string within quotes

   call seek_close_paren( assertif_dir, 1, logical_len )

!  verify found logical expression

   no_logical_exp: if( logical_len == 0 )then

      call msg_quit( "can't find assertif expression: " // trim( assertif_dir))

   endif no_logical_exp

! ----------------------------------------------------------------------

!  if active block, process assertif directive

   active_line: if( if_construct% now_selected )then

!  evaluate logical expression

      call eval_log_expr( assertif_dir( 2: logical_len - 1) // blank, logical_value)

!  if true, process assert

      expr_true: if( logical_value )then

         call process_assert_directive( assertif_dir( logical_len + 1: ))

      endif expr_true

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed assertif directive")

      endif extensions

! ----------------------------------------------------------------------

   endif active_line

! ----------------------------------------------------------------------

!  process_assertif_directive() exit

return

! **********************************************************************

!  process_assertif_directive()

end subroutine process_assertif_directive

! **********************************************************************
! **********************************************************************

!  process_dump_directive() process an dump directive

subroutine process_dump_directive

! **********************************************************************

!  entry: dump_dir is quoted dump condition

!  exit: dumpion code is written to the output

! **********************************************************************

!  process_dump_directive() constants

! ----------------------------------------------------------------------

!  dump banner

   character( len= *), parameter :: dump_header = 'dump coco symbols: file: '

!  mark integers and logicals as constant or variable

   character( len= *), parameter :: variable_str = 'variable'

   character( len= *), parameter :: constant_str = 'constant'

! ----------------------------------------------------------------------

!  column widths

   integer, parameter :: vc_len = max( len( variable_str), len( constant_str) ) + 1

   integer, parameter :: pad_len = 2

! ----------------------------------------------------------------------

!  names of types

   integer, parameter :: type_len = 10

   character( len= type_len), dimension( 5), parameter :: type_str = (/ &
                                                            'integer   ', &
                                                            'logical   ', &
                                                            'macro     ', &
                                                            'text      ', &
                                                            '<unknown> ' /)

   integer, parameter :: unknown_idx = 5

   character( len= *), parameter :: undefined_str = '<undefined>'

! **********************************************************************

!  process_dump_directive() local

! ----------------------------------------------------------------------

!  convert line number to string

   character( len= conversion_len) :: line_number_str

!  point to symbols

   type( symbol_t), pointer :: symbol_ptr

!  point to next character to be added to report line

   integer :: next_char

   integer :: arg_len

   integer :: this_arg

! **********************************************************************

!  process_dump_directive() text

continue

! ----------------------------------------------------------------------

!  if active block, process dump directive

   active_line: if( if_construct% now_selected )then

      write( line_number_str, fmt= conversion_fmt) current_file% lines_transfered

      log_line = dump_header // trim( current_file% name_str) &
                             // ", line: " // trim( adjustl( line_number_str))

      call write_source_line( log_file)

! ----------------------------------------------------------------------

!  write a line for each symbol on the list

      symbol_ptr => first_symbol

      all_symbols: do while( associated( symbol_ptr) )

         log_line = blank

!  line varies by type fo symbol

         which_type: select case( symbol_ptr% type_code)

! ----------------------------------------------------------------------

!  type integer

         case( type_integer) which_type

!  type and name

            log_line( : type_len) = type_str( symbol_ptr% type_code)

            log_line( type_len: type_len + symbol_name_len) = symbol_ptr% name_str

!  constant or variable

            next_char = type_len + symbol_name_len

            integer_v_or_c: if( symbol_ptr% constant )then

               log_line( next_char: next_char + vc_len) = constant_str

            else integer_v_or_c

               log_line( next_char: next_char + vc_len) = variable_str

            endif integer_v_or_c

!  defined or undefined

            next_char = next_char + vc_len + pad_len

!  if defined write the value

            integer_defined: if( symbol_ptr% defined )then

               write( unit= log_line( next_char: next_char + conversion_len), fmt= conversion_fmt) &
                      symbol_ptr% integer_value

               log_line( next_char: next_char + conversion_len + pad_len) = &
                         adjustr( log_line( next_char: next_char + conversion_len + pad_len) )

            else integer_defined

               log_line( next_char + pad_len: ) = undefined_str

            endif integer_defined

! ----------------------------------------------------------------------

!  type logical

         case( type_logical) which_type

!  type and name

            log_line( : type_len) = type_str( symbol_ptr% type_code)

            log_line( type_len: type_len + symbol_name_len) = symbol_ptr% name_str

!  constant or variable

            next_char = type_len + symbol_name_len

            logical_v_or_c: if( symbol_ptr% constant )then

               log_line( next_char: next_char + vc_len) = constant_str

            else logical_v_or_c

               log_line( next_char: next_char + vc_len) = variable_str

            endif logical_v_or_c

!  defined or undefined

            next_char = next_char + vc_len + pad_len

            logical_defined: if( symbol_ptr% defined )then

!  if defined write the value

               logical_t_f: if( symbol_ptr% logical_value )then

                  log_line( next_char: next_char + conversion_len) = true_str

               else logical_t_f

                  log_line( next_char: next_char + conversion_len) = false_str

               endif logical_t_f

               log_line( next_char: next_char + conversion_len + pad_len) = &
                         adjustr( log_line( next_char: next_char + conversion_len + pad_len) )

            else logical_defined

               log_line( next_char + pad_len: ) = undefined_str

            endif logical_defined

! ----------------------------------------------------------------------

!  type macro

         case( type_macro) which_type

!  type and name

            log_line( : type_len) = type_str( symbol_ptr% type_code)

            log_line( type_len: type_len + symbol_name_len) = symbol_ptr% name_str

!  has arg list or not

            next_char = type_len + symbol_name_len

            log_line( next_char: next_char) = open_paren

            next_char = next_char + len( open_paren)

            macro_args: if( associated( symbol_ptr% dummy_args) )then

               arg_len = len_trim( symbol_ptr% dummy_args( 1) )

               log_line( next_char: next_char + arg_len) = blank // symbol_ptr% dummy_args( 1)

               next_char = next_char + arg_len + len( blank)

               next_macro_arg: do this_arg = 2, size( symbol_ptr% dummy_args)

                  arg_len = len_trim( symbol_ptr% dummy_args( this_arg) )

                  log_line( next_char: next_char + arg_len + 1) = &
                          comma // blank // symbol_ptr% dummy_args( this_arg)

                  next_char = next_char + arg_len + len( comma // blank)

               enddo next_macro_arg

            endif macro_args

            log_line( next_char: next_char) = close_paren

!  value

            next_char = max( next_char + pad_len, type_len + symbol_name_len + vc_len + pad_len)

            log_line( next_char: ) = symbol_ptr% macro_value

! ----------------------------------------------------------------------

! type text

         case( type_text) which_type

!  type and name

            log_line( : type_len) = type_str( symbol_ptr% type_code)

            log_line( type_len: type_len + symbol_name_len) = symbol_ptr% name_str

!  has arg list or not

            next_char = type_len + symbol_name_len

            log_line( next_char: next_char) = open_paren

            next_char = next_char + len( open_paren)

            text_args: if( associated( symbol_ptr% dummy_args) )then

               arg_len = len_trim( symbol_ptr% dummy_args( 1) )

               log_line( next_char: next_char + arg_len) = blank // symbol_ptr% dummy_args( 1)

               next_char = next_char + arg_len + len( blank)

               next_text_arg: do this_arg = 2, size( symbol_ptr% dummy_args)

                  arg_len = len_trim( symbol_ptr% dummy_args( this_arg) )

                  log_line( next_char: next_char + arg_len + 1) = comma // blank // symbol_ptr% dummy_args( this_arg)

                  next_char = next_char + arg_len + len( comma // blank)

               enddo next_text_arg

            endif text_args

            log_line( next_char: next_char) = close_paren

!  first line of text

            next_char = max( next_char + pad_len, type_len + symbol_name_len + vc_len + pad_len)

            got_text_lines: if( associated( symbol_ptr% text_lines) )then

               log_line( next_char: ) = symbol_ptr% text_lines( 1)

            endif got_text_lines

! ----------------------------------------------------------------------

!  unknown type

         case default which_type

!  type and name

            log_line( : type_len) = type_str( unknown_idx)

            log_line( type_len: type_len + symbol_name_len) = symbol_ptr% name_str

         end select which_type

!  next_char = next_char + symbol_name_len

         call write_source_line( log_file)

         symbol_ptr => symbol_ptr% next

      enddo all_symbols

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed dump directive")

      endif extensions

   endif active_line

! ----------------------------------------------------------------------

!  process_dump_directive() exit

return

! **********************************************************************

!  process_dump_directive()

end subroutine process_dump_directive

! **********************************************************************
! **********************************************************************

!  process_text_directive() process an text declaration

subroutine process_text_directive( text_dir)

! **********************************************************************

!  process_text_directive() interface

! ----------------------------------------------------------------------

!  a statement containing a text directive

character( len= *), intent( in) :: text_dir                         

! **********************************************************************

!  entry: text_dir is quoted text condition

!  exit: text code is stored in the text variable on the symbol list

! **********************************************************************

!  process_text_directive() local

! ----------------------------------------------------------------------

!  name of symbol

   character( len= symbol_name_len) :: text_name

!  results of decoding statement

   integer :: symbol_len

!  point to next character to be decoded

   integer :: next_char

! **********************************************************************

!  process_text_directive() text

continue

   next_char = 1

! ----------------------------------------------------------------------

!  if active block, process text declaration

   active_line: if( if_construct% now_selected )then

!  extract text name

      call get_text_name( text_dir( next_char: ), text_name, symbol_len)

!  store text in symbol list

      next_char = next_char + symbol_len

      call add_text( text_dir( next_char: ), text_name)

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed text declaration " // trim( text_name) )

      endif extensions

!  count text blocks

      total% text_blocks = total% text_blocks + 1

   endif active_line

! ----------------------------------------------------------------------

!  process_text_directive() exit

return

! **********************************************************************

!  process_text_directive()

end subroutine process_text_directive

! **********************************************************************
! **********************************************************************

!  add_text() copy text block to symbol table

subroutine add_text( text_decl_str, text_name)

! **********************************************************************

!  add_text() interface

! ----------------------------------------------------------------------

!  a statement containing a text declaration

character( len= *), intent( in) :: text_decl_str

!  the name of the text

character( len= *), intent( in) :: text_name                        

! **********************************************************************

!  entry: mac_decl_str is blank_compress_lower_case logical declaration statement past the double colon
!         "" | "=..."
!         text_name is the symbol name
!         is_const is true if this is a constant declaration
!         processing_set_file is true if this was found in the setfile

!  exit: logical declaration is added to the logical symbol list or error exit

! **********************************************************************

!  add_text() constants

! ----------------------------------------------------------------------

!  end of a text block

   character( len= *), parameter :: endtext_str = 'endtext'

! **********************************************************************

!  add_text() local

! ----------------------------------------------------------------------

!  pointer to pre-existing symbol

   type( symbol_t), pointer :: symbol_ptr

!  location of close parenthesis

   integer :: close_idx

!  buffer for text argument list

   character( len= buffer_len) :: arg_list_str

! ----------------------------------------------------------------------

!  the text scratch file

   type( file_t) :: text_file

! ----------------------------------------------------------------------

!  copy buffer

   character( len= buffer_len) :: statement

!  lines of text in the text block

   integer :: this_line

   integer :: this_arg

   character( len= target_len) :: search_str

   integer :: arg_idx

   integer :: arg_len

   logical :: complete

   integer :: next_char

!  allocation status

   integer :: astat

! **********************************************************************

!  add_text() text

continue

! ----------------------------------------------------------------------

!  check if already on text list

   call seek_symbol_name( text_name, symbol_ptr)

   have_a_name: if( associated( symbol_ptr) )then

      duplicate_name: select case( symbol_ptr% type_code)
   
      case( type_text) duplicate_name
   
         call msg_quit( "duplicate text declaration: " // trim( text_name) )

      case( type_integer) duplicate_name
   
         call msg_quit( "text name already defined as an integer: " // trim( text_name) )

      case( type_logical) duplicate_name
   
         call msg_quit( "text name already defined as a logical: " // trim( text_name) )

      case( type_macro) duplicate_name
   
         call msg_quit( "text name already defined as a macro: " // trim( text_name) )

      case default duplicate_name

         call msg_quit( "duplicate name: " // trim( text_name) )

      end select duplicate_name

   endif have_a_name

! ----------------------------------------------------------------------

!  allocate next text

   call add_symbol( text_name, symbol_ptr)

! ----------------------------------------------------------------------

!  build complete text entry

   symbol_ptr% type_code = type_text

!  determine if declaration specifies an arg list

   got_arglist: if( text_decl_str( : len( open_paren)) == open_paren )then

      call seek_close_paren( text_decl_str, len( open_paren), close_idx)

      no_close_paren: if( close_idx == 0 )then

         call msg_quit( "no closing parenthesis: " // trim( text_decl_str) )

      endif no_close_paren

      arg_list_str = text_decl_str( len( open_paren) + 1: close_idx - 1)

      call process_dummy_arglist( arg_list_str, symbol_ptr)

      call verify_dummy_args( symbol_ptr)

   else got_arglist

      nullify( symbol_ptr% dummy_args)

   endif got_arglist

! ----------------------------------------------------------------------

!  store text value (read into scratch file, count lines, allocate storage, copy to storage)

! ----------------------------------------------------------------------

!  initialize the text file variable

   text_file = file_t( text_unit, null_string, &
                             null_string, null(), &
                             0, 0, &
                             .false., .true.)

!  open the set text file

   call open_scratch( text_file)

!  start as if with a complete statement

   complete = .true.

! ----------------------------------------------------------------------

!  main read setfile lines loop

   read_lines: do

! ----------------------------------------------------------------------

!  read a text line from the current source file

      read( unit= current_file% logical_unit, fmt= current_file% format_str, &
            iostat= current_file% io_status) current_file% line

      read_error: if( current_file% io_status > 0 )then

         call msg_quit( "read text failed: " // trim( current_file% name_str))

      endif read_error

! ----------------------------------------------------------------------

!  read until end of file or complete statement

      read_eof: if( current_file% io_status < 0 )then

         total% input_lines = total% input_lines + current_file% lines_transfered

         call msg_quit( "end of file encountered within text block")

      endif read_eof

!  count lines

      current_file% lines_transfered = current_file% lines_transfered + 1

      total% text_lines = total% text_lines + 1

!  write all lines to the output as coco lines

      call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  process coco lines

      coco_line: if( current_file% line( : len( coco_key)) == coco_key )then

!  count coco lines

         total% coco_lines = total% coco_lines + 1

!  ignore coco comments

         coco_statement: if( is_coco_statement( current_file% line( len( coco_key) + 1: )) )then

!  gather a complete statement

            call gather_coco_statement( current_file% line, statement, complete)

!  if incomplete, go get rest of statement

            got_statement: if( .not. complete )then

               cycle read_lines

            endif got_statement

! ----------------------------------------------------------------------

!  check for the end text statement

            end_text: if( statement == endtext_str // symbol_ptr% name_str &
                     .or. statement == endtext_str )then

               exit read_lines

            endif end_text

!  check for certain directives in the text block

            call verify_text_directive( statement)

         endif coco_statement

!  source lines

      else coco_line

         continuation_error: if( .not. complete )then

            call msg_quit( "source line in continued coco statement")

         endif continuation_error

!  end processing text statements

      endif coco_line

!  write the text line

      write( unit= text_file% logical_unit, iostat= text_file% io_status) text_file% line

      write_text: if( text_file% io_status > 0 )then

         call msg_quit( "write text file failed: " // trim( text_file% line))

      endif write_text

! count text lines

      text_file% lines_transfered = text_file% lines_transfered + 1

!  end main read setfile lines loop

   enddo read_lines

! ----------------------------------------------------------------------

!  check for no lines in text block

   null_text: if( text_file% lines_transfered == 0 )then

      call close_scratch( text_file)

      call msg_quit( "text block has no lines: " // trim( symbol_ptr% name_str) )

   endif null_text

! ----------------------------------------------------------------------

!  allocate array for text

   allocate( symbol_ptr% text_lines( text_file% lines_transfered), stat= astat)

   alloc_error: if( astat > 0 )then

      call msg_quit( "allocate text block failed " // trim( text_name) )

   endif alloc_error

!  count text lines defined

   total% text_lines = total% text_lines + 1


!  rewind text scratch file

   rewind( unit= text_file% logical_unit, iostat= text_file% io_status)

   rewind_text: if( text_file% io_status > 0 )then

      call msg_quit( "rewind text scratch file failed")

   endif rewind_text

!  copy text scratch file to array

   copy: do this_line = 1, size( symbol_ptr% text_lines)

      symbol_ptr% text_lines( this_line) = blank

      read( unit= text_file% logical_unit, iostat= text_file% io_status) symbol_ptr% text_lines( this_line)

      read_text: if( text_file% io_status > 0 )then

         call msg_quit( "read text scratch file failed")

      endif read_text

!  force lines to lower case

      each_char: do next_char = 1, len( symbol_ptr% text_lines( this_line))

         to_lower: select case( symbol_ptr% text_lines( this_line)( next_char: next_char))

         case( 'A': 'Z') to_lower

            symbol_ptr% text_lines( this_line)( next_char: next_char) = &
               achar( iachar( symbol_ptr% text_lines( this_line)( next_char: next_char)) + change_case)

         end select to_lower

      enddo each_char

   enddo copy

!  verify whether each dummy arg appears in the text block somewhere

   has_dummy_args: if( associated( symbol_ptr% dummy_args) )then

      check_arg: do this_arg = lbound( symbol_ptr% dummy_args, dim= 1), ubound( symbol_ptr% dummy_args, dim= 1)

         arg_idx = 0

         search_str = arg_key // trim( symbol_ptr% dummy_args( this_arg) ) // arg_key

         arg_len = len_trim( search_str)

         check_line: do this_line = lbound( symbol_ptr% text_lines, dim= 1), ubound( symbol_ptr% text_lines, dim= 1)

            arg_idx = max( arg_idx, index( symbol_ptr% text_lines( this_line), search_str( : arg_len) ) )

         enddo check_line

         not_found: if( arg_idx == 0 )then

            call msg_quit( "dummy arg " // search_str( : arg_len) &
                           // " not found in text: " // trim( symbol_ptr% name_str) )

         endif not_found

      enddo check_arg

   endif has_dummy_args

!  close text scratch file

   call close_scratch( text_file)

! ----------------------------------------------------------------------

!  add_text() exit

return

! **********************************************************************

!  add_text()

end subroutine add_text

! **********************************************************************
! **********************************************************************

!  verify_text_directive() check that no invalid directives appear in a text block

subroutine verify_text_directive( text_stmt)

! **********************************************************************

!  verify_text_directive() interface

! ----------------------------------------------------------------------

!  a statement from a text block

character( len= *), intent( in) :: text_stmt                       

! **********************************************************************

!  entry: text_stmt is a blank_compress_lower_case coco directive past the coco key
!         which must not contain setfile directives or any of:
!         "include..." | "integer..." | "logical..." |
!         "macro..." | "text..." | "copy..."

!  exit: if it exists, the directive is found and flagged

! **********************************************************************

!  verify_text_directive() local

   logical :: flag

! **********************************************************************

!  verify_text_directive() text

continue

! ----------------------------------------------------------------------

!  which directive?

   flag = .false.

! ----------------------------------------------------------------------

!  stop directive

   which_directive: if( text_stmt( : len( include_str)) == include_str )then

      flag = .true.

! ----------------------------------------------------------------------

!  integer declaration

   elseif( text_stmt( : len( integer_str)) == integer_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  integer constant declaration

   elseif( text_stmt( : len( integer_constant_str)) == integer_constant_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  logical declaration

   elseif( text_stmt( : len( logical_str)) == logical_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  logical constant declaration

   elseif( text_stmt( : len( logical_constant_str)) == logical_constant_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  macro declaration

   elseif( text_stmt( : len( macro_str)) == macro_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  text declaration

   elseif( text_stmt( : len( text_str)) == text_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  copy directive

   elseif( text_stmt( : len( copy_str)) == copy_str )then which_directive

      flag = .true.

! ----------------------------------------------------------------------

!  which directive?

   endif which_directive

! ----------------------------------------------------------------------

!  complain and quit if any directives found which shouldn't be in a text block

   dir_check: if( flag )then

      call msg_quit( "illegal directive in text block: " &
                     // trim( text_stmt) )

   endif dir_check

! ----------------------------------------------------------------------

!  verify_text_directive() exit

return

! **********************************************************************

!  verify_text_directive()

end subroutine verify_text_directive

! **********************************************************************
! **********************************************************************

!  process_copy_directive() process a coco copy directive

subroutine process_copy_directive( copy_dir)

! **********************************************************************

!  process_copy_directive() interface

! ----------------------------------------------------------------------

!  a statement containing a copy directive

character( len= *), intent( in) :: copy_dir

! **********************************************************************

!  entry: copy directive

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_copy_directive() constants

! ----------------------------------------------------------------------

!  mark beginning and end of text

   character( len= *), parameter :: begin_txt = '??! TEXT '

   character( len= *), parameter :: end_txt = '??! END TEXT '

! **********************************************************************

!  process_copy_directive() local

! ----------------------------------------------------------------------

!  use name from directive to find text block pointer

   character( len= symbol_name_len) :: text_name

!  length of text block name

   integer :: name_len

!  find beginning of name

   integer :: name_idx

!  find end of name

   integer :: end_name_idx

!  pointer to text block

   type( symbol_t), pointer :: text_ptr

!  gather a coco statement from the text block

   character( len= buffer_len) :: statement

!  loop through the text block lines

   integer :: this_line

!  find open parenthesis on copy directive

   integer :: open_paren_idx

!  find close parenthesis

   integer :: close_paren_idx

!  number of lines in text block or zero if none

   integer :: text_lines_size

!  if args found

   logical :: process_args

!  communicate with gather_statement()

   logical :: complete

! **********************************************************************

!  process_copy_directive() text

continue

! ----------------------------------------------------------------------

!  check for valid directive

   name_idx = 1

   call get_text_name( copy_dir, text_name, name_len)

   end_name_idx = name_idx + name_len - 1

   call get_text_ptr( copy_dir( name_idx: end_name_idx), text_ptr)

   check_name: if( .not. associated( text_ptr) )then

      call msg_quit( "text name not found: " // trim( copy_dir( name_idx: )) )

   endif check_name

! ----------------------------------------------------------------------

!  if active block, process text declaration

   active_line: if( if_construct% now_selected )then

!  test first character after name

      open_paren_idx = end_name_idx + 1

!  check that if text has dummy args, copy has actual args, and vice versa

      have_args: if( associated( text_ptr% dummy_args) )then

!  text with args

         no_args: if( copy_dir( open_paren_idx: open_paren_idx) /= open_paren )then

            call msg_quit( "no actual arguments for text: " // trim( text_ptr% name_str) )

         endif no_args

         call seek_close_paren( copy_dir, open_paren_idx, close_paren_idx)

         process_args = .true.

      else have_args

!  text without args

         got_args: if( copy_dir( open_paren_idx: open_paren_idx) == open_paren )then

            call msg_quit( "no dummy arguments for text: " // trim( copy_dir))

         endif got_args

         process_args = .false.

!  block has/has not args

      endif have_args

! ----------------------------------------------------------------------

!  mark the beginning of the text

      line = begin_txt // text_ptr% name_str

      call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  get number of lines in text block

      have_lines: if( associated( text_ptr% text_lines) )then

         text_lines_size = size( text_ptr% text_lines)

      else have_lines

         text_lines_size = 0

      endif have_lines

! ----------------------------------------------------------------------

!  loop thru text block lines

      copy_lines: do this_line = 1, text_lines_size

         line = text_ptr% text_lines( this_line)

!  coco lines or source lines

         coco_lines: if( line( : len( coco_key)) == coco_key )then

!  write coco line to the output

            call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  process coco lines, ignore coco comments

            coco_construct: if( is_coco_statement( line( len( coco_key) + 1: )) )then

!  gather a complete coco statement

               call gather_coco_statement( line, statement, complete)

!  if not yet a complete statement, get next line

               incomplete: if( .not. complete )then

                  cycle copy_lines

               endif incomplete

!  process (permitted in a block) directives

               call process_block_directive( statement)

               output_file% line => line

            endif coco_construct

! ----------------------------------------------------------------------

!  source lines

         else coco_lines

!  if args substitute in text line

            go_args: if( process_args )then

               text_ptr% macro_value = line

               call process_actual_arglist( copy_dir( open_paren_idx + 1: close_paren_idx - 1), &
                                            line, text_ptr)

            endif go_args

!  editing source lines is enabled

            edit_source_args: if( options% edit_source )then

!  if ? present, edit source line

               edit_line_args: if( index( line, arg_key) > 0 )then

                  call edit_source_line( line)

               endif edit_line_args

            endif edit_source_args

!  finally, write out the line

            call write_source_line( output_file)

         endif coco_lines

      enddo copy_lines

      total% copied_lines = total% copied_lines + text_lines_size

!  mark the end of the text

      line = end_txt // text_ptr% name_str

      call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed copy directive " // trim( text_ptr% name_str) )

      endif extensions

!  process active lines only

   endif active_line

! ----------------------------------------------------------------------

!  process_copy_directive() exit

return

! **********************************************************************

!  process_copy_directive()

end subroutine process_copy_directive

! **********************************************************************
! **********************************************************************

!  process_copyif_directive() process a coco copyif directive

subroutine process_copyif_directive( copyif_dir)

! **********************************************************************

!  process_copyif_directive() interface

! ----------------------------------------------------------------------

!  a statement containing a copy directive

character( len= *), intent( in) :: copyif_dir

! **********************************************************************

!  entry: copy directive

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_copyif_directive() local

! ----------------------------------------------------------------------

!  communicate with seek_close_paren()

   integer :: logical_len

!  communicate with eval_log_expr()

   logical :: logical_value

! **********************************************************************

!  process_copyif_directive() text

continue

! ----------------------------------------------------------------------

!  check syntax- condition must be string within quotes

   call seek_close_paren( copyif_dir, 1, logical_len )

!  verify found logical expression

   no_logical_exp: if( logical_len == 0 )then

      call msg_quit( "can't find copyif expression: " // trim( copyif_dir))

   endif no_logical_exp

! ----------------------------------------------------------------------

!  if active block, process assertif directive

   active_line: if( if_construct% now_selected )then

!  evaluate logical expression

      call eval_log_expr( copyif_dir( 2: logical_len - 1) // blank, logical_value)

!  if true, process copy

      expr_true: if( logical_value )then

         call process_copy_directive( copyif_dir( logical_len + 1: ))

      endif expr_true

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed copyif directive")

      endif extensions

! ----------------------------------------------------------------------

   endif active_line

! ----------------------------------------------------------------------

!  process_copyif_directive() exit

return

! **********************************************************************

!  process_copyif_directive()

end subroutine process_copyif_directive

! **********************************************************************
! **********************************************************************

!  process_block_directive() process a coco text block directive

subroutine process_block_directive( block_stmt)

! **********************************************************************

!  process_block_directive() interface

! ----------------------------------------------------------------------

!  a statement from a text block

character( len= *), intent( in) :: block_stmt

! **********************************************************************

!  entry: coco_stmt is a blank_compress_lower_case coco directive past the coco key
!         "stop..." | "message..." | "if..." | "elseif..." | "else..." |
!         "endif..." | "assert..." | "name=..."

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_block_directive() local

! ----------------------------------------------------------------------

!  point to location of symbol

   type( symbol_t), pointer :: symbol_ptr

!  pointer to equals

   integer :: eq_idx

!  expression string is after the equals

   integer :: expr_idx

! **********************************************************************

!  process_block_directive() text

continue

! ----------------------------------------------------------------------

!  which directive?

! ----------------------------------------------------------------------

!  detect assignment statements assigning to variables named by keywords

      eq_idx = scan( block_stmt( : symbol_name_len + 1), equals)

      got_equals: if( eq_idx > 0 )then

         call seek_symbol_name( block_stmt( : eq_idx - 1), symbol_ptr)

      endif got_equals

! ----------------------------------------------------------------------

!  which directive?

! ----------------------------------------------------------------------

!  assignment directive

   which_directive: if( associated( symbol_ptr) )then

!  up to the equals must be a declared name

      expr_idx = eq_idx + len( equals)

!  must be an integer or logical variable

      integer_or_logical: select case( symbol_ptr% type_code)
         
      case( type_integer) integer_or_logical

         call process_integer_assignment( block_stmt( expr_idx: ), symbol_ptr)

      case( type_logical) integer_or_logical

         call process_logical_assignment( block_stmt( expr_idx: ), symbol_ptr)

      case default integer_or_logical

         call msg_quit( "variable must be an integer or a logical: " // trim( symbol_ptr% name_str) )
         
      end select integer_or_logical

      nullify( symbol_ptr)

! ----------------------------------------------------------------------

!  stop directive

   elseif( block_stmt( : len( stop_str)) == stop_str )then which_directive

      call process_stop_directive( block_stmt( len( stop_str) + 1: ) )

! ----------------------------------------------------------------------

!  message directive

   elseif( block_stmt( : len( message_str)) == message_str )then which_directive

      call process_message_directive( block_stmt( len( message_str) + 1: ) )

! ----------------------------------------------------------------------

!  if directive

   elseif( block_stmt( : len( if_str)) == if_str )then which_directive

      call process_if_directive( block_stmt( len( if_str) + 1: ) )

! ----------------------------------------------------------------------

!  elseif directive

   elseif( block_stmt( : len( elseif_str)) == elseif_str )then which_directive

      call process_elseif_directive( block_stmt( len( elseif_str) + 1: ) )

! ----------------------------------------------------------------------

!  else directive

   elseif( block_stmt( : len( else_str)) == else_str )then which_directive

      call process_else_directive( block_stmt( len( else_str) + 1: ) )

! ----------------------------------------------------------------------

!  endif directive

   elseif( block_stmt( : len( endif_str)) == endif_str )then which_directive

      call process_endif_directive( block_stmt( len( endif_str) + 1: ) )

! ----------------------------------------------------------------------

!  assert declaration

   elseif( block_stmt( : len( assert_str)) == assert_str )then which_directive

      call process_assert_directive( block_stmt( len( assert_str) + 1: ))

! ----------------------------------------------------------------------

!  cannot process this directive

   else which_directive

         call msg_quit( "unknown block directive: " // trim( block_stmt))

! ----------------------------------------------------------------------

!  which directive?

   endif which_directive

! ----------------------------------------------------------------------

!  process_block_directive() exit

return

! **********************************************************************

!  process_block_directive()

end subroutine process_block_directive

! **********************************************************************
! **********************************************************************

!  process_doc_directive() process a coco doc directive

subroutine process_doc_directive( doc_dir)

! **********************************************************************

!  process_doc_directive() interface

! ----------------------------------------------------------------------

!  a statement containing a copy directive

character( len= *), intent( in) :: doc_dir

! **********************************************************************

!  entry: doc directive

!  exit: the directive is processed or error exit

! **********************************************************************

!  process_doc_directive() constants

! ----------------------------------------------------------------------

!  mark beginning and end of text

   character( len= *), parameter :: enddoc_str = 'enddoc'

! **********************************************************************

!  process_doc_directive() local

! ----------------------------------------------------------------------

!  mark beginning and end of text

   character( len= buffer_len) :: statement

   logical :: complete

! **********************************************************************

!  process_doc_directive() text

continue

! ----------------------------------------------------------------------

!  check for valid directive

   valid_doc: if( doc_dir /= blank )then

      call msg_quit( "invalid characters at end of doc directive: " // doc_dir)

   endif valid_doc

! ----------------------------------------------------------------------

!  if active block, process text declaration

   active_line: if( if_construct% now_selected )then

!  if no doc file open, error

      no_doc: if( .not. doc_file% named_file )then

         call msg_quit( 'no document file when doc directive encountered')

      endif no_doc

!  connect doc file to line buffer

      doc_file% line => line

!  read until end doc found or ( enf or error)

      copy_lines: do

!  read doc lines from input

         read( unit= current_file% logical_unit, fmt= current_file% format_str, &
               iostat= current_file% io_status) current_file% line

         read_error: if( current_file% io_status > 0 )then

            call msg_quit( "read input failed: " // trim( current_file% name_str))

         endif read_error

! ----------------------------------------------------------------------

!  read until end of file or complete statement

         read_eof: if( current_file% io_status < 0 )then

            total% input_lines = total% input_lines + current_file% lines_transfered

            call msg_quit( "end of file encountered within doc ... end doc")

         endif read_eof

!  write line to output as coco line

         call write_coco_line( output_file)

!  count current file lines

         current_file% lines_transfered = current_file% lines_transfered + 1

! ----------------------------------------------------------------------

!  process setfile lines or error if source lines

         coco_line: if( line( : len( coco_key)) == coco_key )then

!  count coco lines

            total% coco_lines = total% coco_lines + 1

!  process setfile lines, ignore coco comments

            coco_statement: if( is_coco_statement( line( len( coco_key) + 1: )) )then

! ----------------------------------------------------------------------

!  read a complete statement line by line

               call gather_coco_statement( line, statement, complete)

!  if not yet a complete statement go get the rest of it

               get_statement: if( .not. complete )then

                  cycle copy_lines

               endif get_statement

!  check statement for end doc directive

               got_end_doc: if( statement == enddoc_str )then

                  exit copy_lines

               endif got_end_doc

!  process setfile lines, ignore coco comments

            endif coco_statement

!  process setfile lines, ignore coco comments

         endif coco_line

!  editing source lines is enabled

         edit_source: if( options% edit_source )then

!  if ? present, edit source line

            edit_line: if( index( line, arg_key) > 0 )then

               call edit_source_line( line)

            endif edit_line

         endif edit_source

!  write line to doc file

         call write_source_line( doc_file)

!  enddo

      enddo copy_lines

!  disconnect doc file from line buffer

      nullify( doc_file% line)

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "processed doc ... end doc directive")

      endif extensions

!  process active lines only

   endif active_line

! ----------------------------------------------------------------------

!  process_doc_directive() exit

return

! **********************************************************************

!  process_doc_directive()

end subroutine process_doc_directive

! **********************************************************************
! **********************************************************************

!  process_output_directive() process include output options

subroutine process_output_directive( output_dir)

! **********************************************************************

!  process_output_directive() interface

! ----------------------------------------------------------------------

!  the output directive from the setfile

character( len= *), intent( in) :: output_dir                      

! **********************************************************************

!  entry: output_opt is a output to be added to the list
!         of directories to be searched for inlcude files

!  exit: output_opt is on the list

! **********************************************************************

!  process_output_directive() local

! ----------------------------------------------------------------------

!  the name of the file to be opened

   character( len= filename_len) :: output_name

!  the length of the quoted string

   integer :: quoted_len

!  the length of the unquoted string

   integer :: unquoted_len

! **********************************************************************

!  process_output_directive() text

continue

! ----------------------------------------------------------------------

!  unquote string on directive

   call unquote_string( output_dir, output_name, unquoted_len, quoted_len)

   no_name: if( unquoted_len == 0 .or. quoted_len == 0 )then

      call msg_quit( "no name found on output directive: " // trim( output_dir) )

   endif no_name

!  verify no extra characters beyond name

   extra_chars: if( output_dir( unquoted_len + 1: ) /= blank )then

      call msg_quit( "extra characters after output file name: " // trim( output_dir))

   endif extra_chars

! ----------------------------------------------------------------------

!  if line is an active line

   active_line: if( if_construct% now_selected )then

!  if the output file has content, copy the setfile to it

   made_output: if( output_file% lines_transfered > 0 )then

!  mark the setfile in the output (whether it is present or not)

         line = mark_set_file

         call write_coco_line( output_file)

! ----------------------------------------------------------------------

!  if processed a setfile

         append_set_file: if( set_file% named_file )then

!  copy setfile file to output

            call copy_set_file

!  if processed setfile

         endif append_set_file

      endif made_output

!  close the current output file

      call close_file( output_file)

!  report to logfile

      summary: if( options% print_summary )then

         call write_report

      endif summary

! ----------------------------------------------------------------------

!  set new name and reset line count

      output_file% name_str = output_name

      output_file% lines_transfered = 0

!  open the new output

      call open_file( output_file)

! ----------------------------------------------------------------------

!  report extension

      extensions: if( options% report_extensions )then

         call msg_continue( "opened new output file: " // trim( output_file% name_str) )

      endif extensions

   endif active_line

! ----------------------------------------------------------------------

!  process_output_directive() exit

return

! **********************************************************************

!  process_output_directive()

end subroutine process_output_directive

! **********************************************************************
! **********************************************************************

!  seek_symbol_name() seek symbol on symbol list

subroutine seek_symbol_name( symbol_str, symbol_ptr)

! **********************************************************************

!  seek_symbol_name() interface

! ----------------------------------------------------------------------

!  the name of the symbol being sought

character( len= *), intent( in) :: symbol_str

!  a pointer to the symbol found

type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  seek_symbol_name() text

continue

! ----------------------------------------------------------------------

!  search symbol list

   nullify( symbol_ptr)

   symbol_ptr => first_symbol

   search_list: do while( associated( symbol_ptr) )

      name_match: if( symbol_str == symbol_ptr% name_str )then

         exit search_list

      endif name_match

      symbol_ptr => symbol_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  seek_symbol_name() exit

return

! **********************************************************************

!  seek_symbol_name()

end subroutine seek_symbol_name

! **********************************************************************
! **********************************************************************

!  get_integer_value() seek symbol on symbol list

subroutine get_integer_value( integer_str, return_value)

! **********************************************************************

!  get_integer_value() interface

! ----------------------------------------------------------------------

!  the name of the integer whose value is sought

character( len= *), intent( in) :: integer_str

!  the value of the integer

integer, intent( out) :: return_value                              

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  get_integer_value() local

! ----------------------------------------------------------------------

!  pointer to search symbol list

   type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  get_integer_value() text

continue

! ----------------------------------------------------------------------

!  search symbol list

   symbol_ptr => first_symbol

   search_list: do while( associated( symbol_ptr) )

      name_match: if( integer_str == symbol_ptr% name_str )then

         right_type: select case( symbol_ptr% type_code)

         case( type_integer) right_type

            value_defined: if( symbol_ptr% defined )then

               return_value = symbol_ptr% integer_value

               return

            else value_defined

               call msg_quit( "integer not defined: " // trim( integer_str) )

            endif value_defined

         case default right_type

            call msg_quit( "need integer, not: " // trim( integer_str) )

         end select right_type

      endif name_match

      symbol_ptr => symbol_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  integer not found

   call msg_quit( "unknown integer: " // trim( integer_str) )

! ----------------------------------------------------------------------

!  get_integer_value() exit

return

! **********************************************************************

!  get_integer_value()

end subroutine get_integer_value

! **********************************************************************
! **********************************************************************

!  get_logical_value() seek symbol on symbol list

subroutine get_logical_value( logical_str, return_value)

! **********************************************************************

!  get_logical_value() interface

! ----------------------------------------------------------------------

!  the name of the logical whose value is sought

character( len= *), intent( in) :: logical_str

!  the value of the logical

logical, intent( out) :: return_value                             

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  get_logical_value() local

! ----------------------------------------------------------------------

!  pointer to search symbol list

   type( symbol_t), pointer :: symbol_ptr

! **********************************************************************

!  get_logical_value() text

continue

! ----------------------------------------------------------------------

!  search symbol list

   symbol_ptr => first_symbol

   search_list: do while( associated( symbol_ptr) )

      name_match: if( logical_str == symbol_ptr% name_str )then

         right_type: select case( symbol_ptr% type_code)

         case( type_logical) right_type

            value_defined: if( symbol_ptr% defined )then

               return_value = symbol_ptr% logical_value

               return

            else value_defined

               call msg_quit( "logical not defined: " // trim( logical_str) )

            endif value_defined

         case default right_type

            call msg_quit( "need logical, not: " // trim( logical_str) )

         end select right_type

      endif name_match

      symbol_ptr => symbol_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  logical not found

   call msg_quit( "unknown logical: " // trim( logical_str) )

! ----------------------------------------------------------------------

!  get_logical_value() exit

return

! **********************************************************************

!  get_logical_value()

end subroutine get_logical_value

! **********************************************************************
! **********************************************************************

!  get_next_integer() seek symbol on symbol list

subroutine get_next_integer( integer_ptr)

! **********************************************************************

!  get_next_integer() interface

! ----------------------------------------------------------------------

!  a pointer to the next integer on the symbol list

type( symbol_t), pointer :: integer_ptr

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  get_next_integer() text

continue

! ----------------------------------------------------------------------

!  start at the symbol list head or continue from previous integer

   start_or_continue: if( associated( integer_ptr) )then

      integer_ptr => integer_ptr% next

   else start_or_continue

      integer_ptr => first_symbol

   endif start_or_continue

!  search symbol list for integers

   search_list: do while( associated( integer_ptr) )

      right_type: select case( integer_ptr% type_code)

      case( type_integer) right_type

         return

      end select right_type

      integer_ptr => integer_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  get_next_integer() exit

return

! **********************************************************************

!  get_next_integer()

end subroutine get_next_integer

! **********************************************************************
! **********************************************************************

!  get_next_logical() seek symbol on symbol list

subroutine get_next_logical( logical_ptr)

! **********************************************************************

!  get_next_logical() interface

! ----------------------------------------------------------------------

!  a pointer to the next logical on the symbol list

type( symbol_t), pointer :: logical_ptr

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  get_next_logical() text

continue

! ----------------------------------------------------------------------

!  start at the symbol list head or continue from previous logical

   start_or_continue: if( associated( logical_ptr) )then

      logical_ptr => logical_ptr% next

   else start_or_continue

      logical_ptr => first_symbol

   endif start_or_continue

!  search symbol list for logicals

   search_list: do while( associated( logical_ptr) )

      right_type: select case( logical_ptr% type_code)

      case( type_logical) right_type

         return

      end select right_type

      logical_ptr => logical_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  get_next_logical() exit

return

! **********************************************************************

!  get_next_logical()

end subroutine get_next_logical                                    

! **********************************************************************
! **********************************************************************

!  get_next_macro() seek symbol on symbol list

subroutine get_next_macro( macro_ptr)

! **********************************************************************

!  get_next_macro() interface

! ----------------------------------------------------------------------

!  a pointer to the next macro on the symbol list

type( symbol_t), pointer :: macro_ptr

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  get_next_macro() text

continue

! ----------------------------------------------------------------------

!  start at the symbol list head or continue from previous macro

   start_or_continue: if( associated( macro_ptr) )then

      macro_ptr => macro_ptr% next

   else start_or_continue

      macro_ptr => first_symbol

   endif start_or_continue

!  search symbol list for macros

   search_list: do while( associated( macro_ptr) )

      right_type: select case( macro_ptr% type_code)

      case( type_macro) right_type

         return

      end select right_type

      macro_ptr => macro_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  get_next_macro() exit

return

! **********************************************************************

!  get_next_macro()

end subroutine get_next_macro

! **********************************************************************
! **********************************************************************

!  get_text_ptr() seek symbol on symbol list

subroutine get_text_ptr( text_str, text_ptr)

! **********************************************************************

!  get_text_ptr() interface

! ----------------------------------------------------------------------

!  the name of the text being sought

character( len= *), intent( in) :: text_str

!  a pointer to the text

type( symbol_t), pointer :: text_ptr                                

! **********************************************************************

!  entry: symbol_str is blank_compress_lower_case logical symbol directive
!         "name..."

!  exit: symbol found or not in logical symbol array

! **********************************************************************

!  get_text_ptr() text

continue

! ----------------------------------------------------------------------

!  start at the symbol list head

   nullify( text_ptr)

   text_ptr => first_symbol

!  search symbol list for texts

   search_list: do while( associated( text_ptr) )

      right_name: if( text_ptr% name_str == text_str )then

         right_type: select case( text_ptr% type_code)

         case( type_text) right_type

            return

         case default right_type

            call msg_quit( "not a text name: " // trim( text_str) )

         end select right_type

      endif right_name

      text_ptr => text_ptr% next

   enddo search_list

! ----------------------------------------------------------------------

!  get_text_ptr() exit

return

! **********************************************************************

!  get_text_ptr()

end subroutine get_text_ptr

! **********************************************************************
! **********************************************************************

!  %%% diagnose and evaluate expressions

! **********************************************************************
! **********************************************************************

!  integer_or_logical() determine type of expression

subroutine integer_or_logical( expr_str, flag)

! **********************************************************************

!  integer_or_logical() interface

! ----------------------------------------------------------------------

!  an expression whose type is to be assertained

character( len= *), intent( in) :: expr_str

!  true if the type is integer

logical, intent( out) :: flag                                       

! **********************************************************************

!  entry: symbol_str is string "..."

!  exit: flag is true if string is an integer expression and false otherwise

! **********************************************************************

!  integer_or_logical() constants

! ----------------------------------------------------------------------

!  search for a character which must be part of a logical expression

   character( len= *), parameter :: logical_chars = '.<>='

!  search for a character which may be part of an integer expression

   character( len= *), parameter :: integer_chars = '+-*/'

! **********************************************************************

!  integer_or_logical() local

! ----------------------------------------------------------------------

!  search results

   integer :: char_idx

   type( symbol_t), pointer ::  symbol_ptr

 ! **********************************************************************

!  integer_or_logical() text

continue

! ----------------------------------------------------------------------

!  does string contain a character which is only in logical expressions?

   char_idx = scan( expr_str, logical_chars)

   got_dot: if( char_idx > 0 )then

      flag = .false.

      return

   endif got_dot

!  does string contain a character which is only in integer expressions?

   char_idx = scan( expr_str, integer_chars)

   got_op: if( char_idx > 0 )then

      flag = .true.

      return

   endif got_op

! ----------------------------------------------------------------------

!  is string an integer or a logical symbol name?

   char_idx = verify( expr_str, alphanum_chars)

   got_name: if( char_idx == 0 )then

      call seek_symbol_name( expr_str, symbol_ptr)

      int_name: if( associated( symbol_ptr) )then

         flag = symbol_ptr% type_code == type_integer

         return

      endif int_name

   endif got_name

! ----------------------------------------------------------------------

!  is string all digits?

   char_idx = verify( expr_str, digit_chars)

   got_digits: if( char_idx == 0 )then

      flag = .true.

      return

   endif got_digits

! ----------------------------------------------------------------------

!  can't classify the expression, punt

   call msg_quit( "can't classify: " // trim( expr_str) )

! ----------------------------------------------------------------------

!  integer_or_logical() exit

return

! **********************************************************************

!  integer_or_logical()

end subroutine integer_or_logical

! **********************************************************************
! **********************************************************************

!  eval_int_expr() evaluate int_expr as an integer

recursive subroutine eval_int_expr( int_expr, value)

! **********************************************************************

!  eval_int_expr() interface

! ----------------------------------------------------------------------

!  the integer expression to be evaluated

character( len= *), intent( in) :: int_expr

!  the value of the integer expression

integer, intent( out) :: value                                      

! **********************************************************************

!  entry: int_expr is blank_compress_lower_case integer int_expr

!  exit: true if value is int_expr value, false otherwise

! **********************************************************************

!  eval_int_expr() constants

! ----------------------------------------------------------------------

!  addition operators

   integer, parameter :: add_op_len = max( len( plus), len( minus) )

!  multiplication operators

   character( len= *), parameter :: times = '*'

   character( len= *), parameter :: divby = '/'

   integer, parameter :: mul_op_len = max( len( times), len( divby) )

!  length of operators

   integer, parameter :: op_len = max( len( plus), len( minus), len( times), len( divby) )

! **********************************************************************

!  eval_int_expr() local

! ----------------------------------------------------------------------

!  operations to be done

   character( len= add_op_len) :: add_op

   character( len= mul_op_len) :: mul_op

! ----------------------------------------------------------------------

!  next operation

   character( len= op_len) :: next_op

! ----------------------------------------------------------------------

!  partial values of the int_expr

   integer :: l_add, r_add

   integer :: l_mul, r_mul

! ----------------------------------------------------------------------

!  pointers to characters

   integer :: next_char

   integer :: next_op_idx

   integer :: expr_len

   integer :: primary_len

! **********************************************************************

!  eval_int_expr() text

continue

! ----------------------------------------------------------------------

!  limits of scan

   next_char = 1

   expr_len = len_trim( int_expr)

!  initialize adds

   add_op = plus

   l_add = 0

! ----------------------------------------------------------------------

!  scan thru int_expr

   add_ops: do while( next_char <= expr_len)

!  find a primary

      call eval_int_primary( int_expr( next_char: ), primary_len, r_add)

      next_op_idx = next_char + primary_len

!  find next operator or end of expression

      add_end: if( next_op_idx <= expr_len )then

         next_op = int_expr( next_op_idx: next_op_idx)

         next_char = next_op_idx + 1

      else add_end

         next_op = blank

         next_char = next_op_idx

      endif add_end

! ----------------------------------------------------------------------

!  initialize for a set of mul ops

      mul_op = next_op

      l_mul = r_add

! ----------------------------------------------------------------------

!  process a set of mul ops

      mul_ops: do while( next_op == times .or. next_op == divby)

!  find a primary

         call eval_int_primary( int_expr( next_char: ), primary_len, r_mul)

         next_op_idx = next_char + primary_len

!  find next operator or end of expression

         mul_end: if( next_op_idx <= expr_len )then

            next_op = int_expr( next_op_idx: next_op_idx)

            next_char = next_op_idx + 1

         else mul_end

            next_op = blank

            next_char = next_op_idx

         endif mul_end

!  do the pending add op

         mul_div: select case( mul_op)

         case( times) mul_div

            l_mul = l_mul * r_mul

         case( divby) mul_div

            l_mul = l_mul / r_mul

         end select mul_div

         mul_op = next_op

      enddo mul_ops

!  product is the right operand

      r_add = l_mul

! ----------------------------------------------------------------------

!  do the pending add op

      add_sub: select case( add_op)

      case( blank, plus) add_sub

         l_add = l_add + r_add

      case( minus) add_sub

         l_add = l_add - r_add

      case default add_sub

         call msg_quit( "unknown arithmetic operator: " // add_op)

      end select add_sub

      add_op = next_op

   enddo add_ops

! ----------------------------------------------------------------------

!  value of integer expression

   value = l_add

! ----------------------------------------------------------------------

!  eval_int_expr() exit

return

! **********************************************************************

!  eval_int_expr()

end subroutine eval_int_expr

! **********************************************************************
! **********************************************************************

!  eval_log_expr() expression is evaluated as a logical

recursive subroutine eval_log_expr( log_expr, value)

! **********************************************************************

!  eval_log_expr() interface

! ----------------------------------------------------------------------

!  the logical expression to be evaluated

character( len= *), intent( in) :: log_expr

!  the value of the expression

logical, intent( out) :: value                                      

! **********************************************************************

!  entry: expression is blank_compress_lower_case logical expression

!  exit: value is expression value

! **********************************************************************

!  eval_log_expr() constants

   integer, parameter :: eqv_op_len = max( len( eqv_str), len( neqv_str))

!  length of the next operator

   integer, parameter :: next_op_len = max( len( or_str), len( and_str), len( eqv_str), len( neqv_str))

! **********************************************************************

!  eval_log_expr() local

! ----------------------------------------------------------------------

!  the current eqv operator

   character( len= eqv_op_len) :: eqv_op

!  the next operator

   character( len= next_op_len) :: next_op

! ----------------------------------------------------------------------

!  point to characters not yet decoded

   integer :: next_char

   integer :: next_op_idx

   integer :: expr_len

   integer :: primary_len

!  false if and but no eqv

   logical :: do_or

!  expression values

   logical :: l_eqv, l_or, l_and

   logical :: r_eqv, r_or, r_and

! **********************************************************************

!  eval_log_expr() text

continue

! ----------------------------------------------------------------------

!  limits of scan

   next_char = 1

   expr_len = len_trim( log_expr)

!  initialize equivalences

   eqv_op = eqv_str

   l_eqv = .true.

! ----------------------------------------------------------------------

!  scan thru log_expr

   eqv_ops: do while( next_char <= expr_len)

!  find a primary and return its length and value

      call eval_log_primary( log_expr( next_char: ), primary_len, r_eqv)

      next_op_idx = next_char + primary_len

!  find next operator or end of expression

      eqv_or_end: if( next_op_idx <= expr_len )then

!  decode which operator

         eqv_next_op: if( log_expr( next_op_idx: next_op_idx + len( eqv_str) - 1) == eqv_str )then

            next_op = log_expr( next_op_idx: next_op_idx + len( eqv_str) - 1)

            next_char = next_op_idx + len( eqv_str)

         elseif( log_expr( next_op_idx: next_op_idx + len( neqv_str) - 1) == neqv_str )then eqv_next_op

            next_op = log_expr( next_op_idx: next_op_idx + len( neqv_str) - 1)

            next_char = next_op_idx + len( neqv_str)

         elseif( log_expr( next_op_idx: next_op_idx + len( or_str) - 1) == or_str )then eqv_next_op

            next_op = log_expr( next_op_idx: next_op_idx + len( or_str) - 1)

            next_char = next_op_idx + len( or_str)

         elseif( log_expr( next_op_idx: next_op_idx + len( and_str) - 1) == and_str )then eqv_next_op

            next_op = log_expr( next_op_idx: next_op_idx + len( and_str) - 1)

            next_char = next_op_idx + len( and_str)

         else eqv_next_op

            call msg_quit( "unknown logical operator: " // trim( log_expr( next_op_idx: ) ))

         endif eqv_next_op

!  past end of expression

      else eqv_or_end

         next_op = blank

         next_char = next_op_idx

      endif eqv_or_end

! ----------------------------------------------------------------------

!  initialize for a set of or ops

      l_or = r_eqv

! ----------------------------------------------------------------------

!  process a set of and ops

      or_ops: do while( next_op == or_str .or. next_op == and_str)

         do_or = next_op == or_str

         or_next: select case( do_or)

         case( .true.) or_next

!  find a primary and return its length and value

            call eval_log_primary( log_expr( next_char: ), primary_len, r_or)

            next_op_idx = next_char + primary_len

!  find next operator or end of expression

            or_end: if( next_op_idx <= expr_len )then

!  decode which operator

               or_next_op: if( log_expr( next_op_idx: next_op_idx + len( eqv_str) - 1) == eqv_str )then

                  next_op = log_expr( next_op_idx: next_op_idx + len( eqv_str) - 1)

                  next_char = next_op_idx + len( eqv_str)

               elseif( log_expr( next_op_idx: next_op_idx + len( neqv_str) - 1) == neqv_str )then or_next_op

                  next_op = log_expr( next_op_idx: next_op_idx + len( neqv_str) - 1)

                  next_char = next_op_idx + len( neqv_str)

               elseif( log_expr( next_op_idx: next_op_idx + len( or_str) - 1) == or_str )then or_next_op

                  next_op = log_expr( next_op_idx: next_op_idx + len( or_str) - 1)

                  next_char = next_op_idx + len( or_str)

               elseif( log_expr( next_op_idx: next_op_idx + len( and_str) - 1) == and_str )then or_next_op

                  next_op = log_expr( next_op_idx: next_op_idx + len( and_str) - 1)

                  next_char = next_op_idx + len( and_str)

               else or_next_op

                  call msg_quit( "unknown logical operator: " // trim( log_expr( next_op_idx: ) ) )

               endif or_next_op

            else or_end

               next_op = blank

               next_char = next_op_idx

            endif or_end

         case default or_next

            r_or = l_or

         end select or_next

! ----------------------------------------------------------------------

!  initialize for a set of and ops

         l_and = r_or

! ----------------------------------------------------------------------

!  process a set of and ops

         and_ops: do while( next_op == and_str)

!  find a primary

            call eval_log_primary( log_expr( next_char: ), primary_len, r_and)

            next_op_idx = next_char + primary_len

!  find next operator or end of expression

            and_end: if( next_op_idx <= expr_len )then

!  decode which operator

               and_next_op: if( log_expr( next_op_idx: next_op_idx + len( eqv_str) - 1) == eqv_str )then

                  next_op = log_expr( next_op_idx: next_op_idx + len( eqv_str) - 1)

                  next_char = next_op_idx + len( eqv_str)

               elseif( log_expr( next_op_idx: next_op_idx + len( neqv_str) - 1) == neqv_str )then and_next_op

                  next_op = log_expr( next_op_idx: next_op_idx + len( neqv_str) - 1)

                  next_char = next_op_idx + len( neqv_str)

               elseif( log_expr( next_op_idx: next_op_idx + len( and_str) - 1) == and_str )then and_next_op

                  next_op = log_expr( next_op_idx: next_op_idx + len( and_str) - 1)

                  next_char = next_op_idx + len( and_str)

               elseif( log_expr( next_op_idx: next_op_idx + len( or_str) - 1) == or_str )then and_next_op

                  next_op = log_expr( next_op_idx: next_op_idx + len( or_str) - 1)

                  next_char = next_op_idx + len( or_str)

               else and_next_op

                  call msg_quit( "unknown logical operator: " // trim( log_expr( next_op_idx: ) ) )

               endif and_next_op

            else and_end

               next_op = blank

               next_char = next_op_idx

            endif and_end

!  do the pending and op

            l_and = l_and .and. r_and

         enddo and_ops

!  product is the right operand

         r_or = l_and

! ----------------------------------------------------------------------

!  do the pending or op

         this_or: select case( do_or)

         case( .true.) this_or

            l_or = l_or .or. r_or

         case default this_or

            l_or = r_or

         end select this_or

      enddo or_ops

!  product is the right operand

      r_eqv = l_or

! ----------------------------------------------------------------------

!  do the pending eqv op

      eqv_neqv: select case( eqv_op)

      case( blank, eqv_str) eqv_neqv

         l_eqv = l_eqv .eqv. r_eqv

      case( neqv_str) eqv_neqv

         l_eqv = l_eqv .neqv. r_eqv

      end select eqv_neqv

      eqv_op = next_op

   enddo eqv_ops

! ----------------------------------------------------------------------

   value = l_eqv

! ----------------------------------------------------------------------

!  eval_log_expr() exit

return

! **********************************************************************

!  eval_log_expr()

end subroutine eval_log_expr

! **********************************************************************
! **********************************************************************

!  eval_rel_expr() a relational expression is evaluated as a logical

subroutine eval_rel_expr( rel_expr, value)

! **********************************************************************

!  eval_rel_expr() interface

! ----------------------------------------------------------------------

!  the relational expression ot be evaluated

character( len= *), intent( in) :: rel_expr

!  the value of the relational expression

logical, intent( out) :: value                                      

! **********************************************************************

!  entry: expression is blank_compress_lower_case relational expression

!  exit: value is expression value

! **********************************************************************

!  eval_rel_expr() local

! ----------------------------------------------------------------------

!  index of symbol entry

   integer :: dot_idx

   integer :: eq_idx, ne_idx, gt_idx, ge_idx, le_idx, lt_idx

   integer :: l_val, r_val

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  eval_rel_expr() text

continue

! ----------------------------------------------------------------------

   dot_idx = index( rel_expr, dot)

! ----------------------------------------------------------------------

!  find a dot?

   got_dot: if( dot_idx > 0 )then

!  seek all operators with dot

      eq_idx = index( rel_expr, dot_eq)

      ne_idx = index( rel_expr, dot_ne)

      gt_idx = index( rel_expr, dot_gt)

      ge_idx = index( rel_expr, dot_ge)

      le_idx = index( rel_expr, dot_le)

      lt_idx = index( rel_expr, dot_lt)

! ----------------------------------------------------------------------

!  find one

      dot_rel_op: if( eq_idx > 0 )then

         expr_str = rel_expr( : eq_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( eq_idx + len( dot_eq): )

         call eval_int_expr( expr_str, r_val)

         value = l_val == r_val

      elseif( ne_idx > 0 )then dot_rel_op

         expr_str = rel_expr( : ne_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( ne_idx + len( dot_ne): )

         call eval_int_expr( expr_str, r_val)

         value = l_val /= r_val

      elseif( ge_idx > 0 )then dot_rel_op

         expr_str = rel_expr( : ge_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( ge_idx + len( dot_ge): )

         call eval_int_expr( expr_str, r_val)

         value = l_val >= r_val

      elseif( le_idx > 0 )then dot_rel_op

         expr_str = rel_expr( : le_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( le_idx + len( dot_le): )

         call eval_int_expr( expr_str, r_val)

         value = l_val <= r_val

      elseif( gt_idx > 0 )then dot_rel_op

         expr_str = rel_expr( : gt_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( gt_idx + len( dot_gt): )

         call eval_int_expr( expr_str, r_val)

         value = l_val > r_val

      elseif( lt_idx > 0 )then dot_rel_op

         expr_str = rel_expr( : lt_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( lt_idx + len( dot_lt): )

         call eval_int_expr( expr_str, r_val)

         value = l_val < r_val

! ----------------------------------------------------------------------

!  unknown relational operator

      else dot_rel_op

         call msg_quit( "no relational operator (.eq., .ne., .gt., .ge., .le., .lt.): " // rel_expr)

      endif dot_rel_op

! ----------------------------------------------------------------------

!  operator without dot

   else got_dot

!  seek all comparison ops

      eq_idx = index( rel_expr, ch_eq)

      ne_idx = index( rel_expr, ch_ne)

      gt_idx = index( rel_expr, ch_gt)

      ge_idx = index( rel_expr, ch_ge)

      le_idx = index( rel_expr, ch_le)

      lt_idx = index( rel_expr, ch_lt)

! ----------------------------------------------------------------------

!  find one

      ch_rel_op: if( eq_idx > 0 )then

         expr_str = rel_expr( : eq_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( eq_idx + len( ch_eq): )

         call eval_int_expr( expr_str, r_val)

         value = l_val == r_val

      elseif( ne_idx > 0 )then ch_rel_op

         expr_str = rel_expr( : ne_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( ne_idx + len( ch_ne): )

         call eval_int_expr( expr_str, r_val)

         value = l_val /= r_val

      elseif( ge_idx > 0 )then ch_rel_op

         expr_str = rel_expr( : ge_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( ge_idx + len( ch_ge): )

         call eval_int_expr( expr_str, r_val)

         value = l_val >= r_val

      elseif( le_idx > 0 )then ch_rel_op

         expr_str = rel_expr( : le_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( le_idx + len( ch_le): )

         call eval_int_expr( expr_str, r_val)

         value = l_val <= r_val

      elseif( gt_idx > 0 )then ch_rel_op

         expr_str = rel_expr( : gt_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( gt_idx + len( ch_gt): )

         call eval_int_expr( expr_str, r_val)

         value = l_val > r_val

      elseif( lt_idx > 0 )then ch_rel_op

         expr_str = rel_expr( : lt_idx - 1)

         call eval_int_expr( expr_str, l_val)

         expr_str = rel_expr( lt_idx + len( ch_lt): )

         call eval_int_expr( expr_str, r_val)

         value = l_val < r_val

! ----------------------------------------------------------------------

!  unknown relational operator

      else ch_rel_op

         call msg_quit( "no relational operator (==, /=, >, >=, <=, <): " // rel_expr)

      endif ch_rel_op

   endif got_dot

! ----------------------------------------------------------------------

!  eval_rel_expr() exit

return

! **********************************************************************

!  eval_rel_expr()

end subroutine eval_rel_expr

! **********************************************************************
! **********************************************************************

!  seek_log_primary() a relational expression is evaluated as a logical

subroutine seek_log_primary( log_expr, op_idx, rel_op_idx)

! **********************************************************************

!  seek_log_primary() interface

! ----------------------------------------------------------------------

!  the logical primary to be evaluated

character( len= *), intent( in) :: log_expr

!  the index of the next operator or end of line after the primary

integer, intent( out) :: op_idx

!  the index of the next relational operator or zero

integer, intent( out) :: rel_op_idx                               

! **********************************************************************

!  entry: find log op before first (if any) open paren or after matching

!  exit: length to first log op

! **********************************************************************

!  seek_log_primary() local

integer :: paren_level

! **********************************************************************

!  seek_log_primary() text

continue

!  initialize while loop parameters

   op_idx = 1

   paren_level = 0

   rel_op_idx = 0

! ----------------------------------------------------------------------

!  scan through expression

   scan_stmt: do while( op_idx <=  len_trim( log_expr))

!  check each character

      which_char: select case( log_expr( op_idx: op_idx))

!  need to track parenthesis level

      case( open_paren) which_char

         paren_level = paren_level + 1

         op_idx = op_idx + len( open_paren)

         cycle scan_stmt

      case( close_paren) which_char

         paren_level = paren_level - 1

         op_idx = op_idx + len( close_paren)

         cycle scan_stmt

      case( dot) which_char

         log_op_at_level_zero: if( paren_level == 0 )then

!  find logical operator

            find_log_op: if( log_expr( op_idx: op_idx + len( or_str) - 1) == or_str )then

               exit scan_stmt

            elseif( log_expr( op_idx: op_idx + len( and_str) - 1) == and_str )then find_log_op

               exit scan_stmt

            elseif( log_expr( op_idx: op_idx + len( eqv_str) - 1) == eqv_str )then find_log_op

               exit scan_stmt

            elseif( log_expr( op_idx: op_idx + len( neqv_str) - 1) == neqv_str )then find_log_op

               exit scan_stmt

            endif find_log_op

         endif log_op_at_level_zero

      end select which_char

!  check for relational operator (which diagnoses a relational expression)

      rel_op_at_level_zero: if( paren_level == 0 )then

         found_rel_op: if( log_expr( op_idx: op_idx + len( dot_eq) - 1) == dot_eq &
                      .or. log_expr( op_idx: op_idx + len( dot_ne) - 1) == dot_ne &
                      .or. log_expr( op_idx: op_idx + len( dot_lt) - 1) == dot_lt &
                      .or. log_expr( op_idx: op_idx + len( dot_le) - 1) == dot_le &
                      .or. log_expr( op_idx: op_idx + len( dot_ge) - 1) == dot_ge &
                      .or. log_expr( op_idx: op_idx + len( dot_gt) - 1) == dot_gt &
                      .or. log_expr( op_idx: op_idx + len( ch_eq) - 1) == ch_eq &
                      .or. log_expr( op_idx: op_idx + len( ch_ne) - 1) == ch_ne &
                      .or. log_expr( op_idx: op_idx + len( ch_lt) - 1) == ch_lt &
                      .or. log_expr( op_idx: op_idx + len( ch_le) - 1) == ch_le &
                      .or. log_expr( op_idx: op_idx + len( ch_ge) - 1) == ch_ge &
                      .or. log_expr( op_idx: op_idx + len( ch_gt) - 1) == ch_gt )then

            rel_op_idx = op_idx

         endif found_rel_op

      endif rel_op_at_level_zero

!  catch unbalanced parenthesis in logical expression

      unbalanced_parens: if( paren_level < 0 )then

         call msg_quit( "unbalanced parenthesis in expression: " // trim( log_expr) )

      endif unbalanced_parens

!  scan next character

      op_idx = op_idx + 1

   enddo scan_stmt

!  point to last character in primary

   op_idx = op_idx - 1

! ----------------------------------------------------------------------

!  seek_log_primary() exit

return

! **********************************************************************

!  seek_log_primary()

end subroutine seek_log_primary

! **********************************************************************
! **********************************************************************

!  eval_int_primary() decode a string to get an integer value

recursive subroutine eval_int_primary( primary_str, primary_len, value)

! **********************************************************************

!  eval_int_primary() interface

! ----------------------------------------------------------------------

!  the integer primary to be evaluated

character( len= *), intent( in) :: primary_str

!  the length of the inetger primary

integer, intent( out) :: primary_len

!  the value of the primary

integer, intent( out) :: value

! **********************************************************************

!  entry: primary_str is a string containing a literal integer

!  exit: primary_len is the length decoded, value is integer value

! **********************************************************************

!  eval_int_primary() local

! ----------------------------------------------------------------------

!  decode integer literal via internal read

   integer :: conversion_iostat

!  process sign separately

   integer :: isign

!  pointers to characters

   integer :: next_char

   integer :: char_idx

   integer :: match_paren

!  decode digit strings

   character( len= conversion_len) :: conversion_str

!  string containing expressions

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  eval_int_primary() text

continue

! ----------------------------------------------------------------------

!  decode unary operator

   next_char = 1

!  evaluate the primary using the expression string

   expr_str = primary_str

! ----------------------------------------------------------------------

!  test first character is minus

   process_sign: select case( expr_str( next_char: next_char) )

! ----------------------------------------------------------------------

   case( minus) process_sign

      next_char = next_char + len( minus)

      primary_len = len( minus)

      isign = -1

! ----------------------------------------------------------------------

!  test first character is plus

   case( plus) process_sign

      next_char = next_char + len( plus)

      primary_len = len( plus)

      isign = 1

! ----------------------------------------------------------------------

!  test first character is neither plus nor minus

   case default process_sign

      primary_len = 0

      isign = 1

   end select process_sign

! ----------------------------------------------------------------------

!  find the value of a variable, a literal, or a parenthesized primary_str

   get_value: select case( expr_str( next_char: next_char) )

! ----------------------------------------------------------------------

!  get the value from the variable

   case( 'a': 'z') get_value

!  seek the value of the symbol name

      char_idx = verify( expr_str( next_char: ) // blank, alphanum_chars) + next_char - 2

      call get_integer_value( expr_str( next_char: char_idx), value)

!  processed the alphanumeric characters

      primary_len = primary_len + char_idx

! ----------------------------------------------------------------------

!  get the value of a literal

   case( '0': '9') get_value

!  find the first character which is not a digit

      char_idx = verify( expr_str( next_char: ) // blank, digit_chars) + next_char - 2

!  decode digits

      conversion_str = expr_str( next_char: char_idx)

      conversion_str = adjustr( conversion_str)

      read( unit= conversion_str, fmt= conversion_fmt, iostat= conversion_iostat) value

!  check read error

      decode: if( conversion_iostat > 0 )then

         call msg_quit( "can't decode: " // primary_str)

      endif decode

!  processed the digit string

      primary_len = primary_len + char_idx

! ----------------------------------------------------------------------

!  get the value of an primary_str

   case( open_paren) get_value

      call seek_close_paren( expr_str, next_char, match_paren)

      found_match: if(  match_paren <= len_trim( primary_str) )then

!  go evaluate the nested expression

         expr_str = primary_str( next_char + 1: match_paren - 1)

         call eval_int_expr( expr_str, value)

!  unmatched parenthesis so complain and quit

      else found_match

         call msg_quit( "unmatched parenthesis: " // trim( primary_str))

      endif found_match

!  processed up to the closing parenthesis

      primary_len = match_paren

! ----------------------------------------------------------------------

!  error: cannot get the value

   case default get_value

      call msg_quit( "bad integer expression: " // trim( primary_str) )

   end select get_value

! ----------------------------------------------------------------------

!  apply sign

   value = value * isign

! ----------------------------------------------------------------------

!  eval_int_primary() exit

return

! **********************************************************************

!  eval_int_primary()

end subroutine eval_int_primary                                    

! **********************************************************************
! **********************************************************************

!  eval_log_primary() decode a string to get an logical value

recursive subroutine eval_log_primary( primary_str, primary_len, value)

! **********************************************************************

!  eval_log_primary() interface

! ----------------------------------------------------------------------

!  the logical primary to be evaluated

character( len= *), intent( in) :: primary_str

!  the length of the logical primary

integer, intent( out) :: primary_len

!  the value of the logical primary

logical, intent( out) :: value                                       

! **********************************************************************

!  entry: primary_str is a string containing a literal logical

!  exit: value is logical value

! **********************************************************************

!  eval_log_primary() local

! ----------------------------------------------------------------------

!  logical "sign"

   logical :: lsign

   integer :: rel_op_idx

!  next character to be decoded

   integer :: next_char

   integer :: match_paren

   character( len= buffer_len) :: expr_str

! **********************************************************************

!  eval_log_primary() text

continue

! ----------------------------------------------------------------------

!  find length of primary and whether it is a relational expression

   call seek_log_primary( primary_str, primary_len, rel_op_idx)

!  decode unary operator

   next_char = 1

! ----------------------------------------------------------------------

!  expression too short to contain a .not.

   process_sign: if( primary_len <= len( not_str) )then

      lsign = .true.

!  expression has a .not.

   elseif( primary_str( next_char: len( not_str)) == not_str )then process_sign

      next_char = next_char + len( not_str)

      lsign = .false.

!  no .not.

   else process_sign

      lsign = .true.

   endif process_sign

! ----------------------------------------------------------------------

!  a logical primary is either a logical expression or a relational expression

   log_or_rel: if( rel_op_idx == 0 )then

! ----------------------------------------------------------------------

!  find the value of a variable, a literal, or a parenthesized expression

      get_value: select case( primary_str( next_char: next_char) )

! ----------------------------------------------------------------------

!  get the value from the variable

      case( 'a': 'z') get_value

!  check whether it's a logical name or error

         call get_logical_value( primary_str( next_char: primary_len), value)

! ----------------------------------------------------------------------

!  get the value of a literal

      case( dot) get_value

!  decode literal value

         literal_value: if( primary_str( next_char: next_char + len( true_str) - 1) == true_str )then

!  found a .true. string

            value = .true.

         elseif( primary_str( next_char: next_char + len( false_str) - 1) == false_str )then literal_value

!  found a .false. string

            value = .false.

!  complain and quit

         else literal_value

            call msg_quit( "bad logical literal: " // trim( primary_str) )

         endif literal_value

! ----------------------------------------------------------------------

!  get the value of an expression

      case( open_paren) get_value

!  seek the closing parenthesis

         call seek_close_paren( primary_str, next_char, match_paren)

!  if found, determine whether it is a logical or (part of a) relational expression

         found_match: if( match_paren <= len_trim( primary_str) )then

!  evaluate the logical expression within parenthesis

            expr_str = primary_str( next_char + 1: match_paren - 1)

            call eval_log_expr( expr_str, value)

!  unmatched parenthesis so complain and quit

         else found_match

            call msg_quit( "unmatched parenthesis: " // trim( primary_str))

         endif found_match

! ----------------------------------------------------------------------

!  error: can't decode logical value

      case default

         call msg_quit( "bad logical primary: " // trim( primary_str))

      end select get_value

! ----------------------------------------------------------------------

!  evaluate the relational expression

   else log_or_rel

         call eval_rel_expr( primary_str( next_char: primary_len), value)

   endif log_or_rel

! ----------------------------------------------------------------------

!  apply sign

   value = value .eqv. lsign

! ----------------------------------------------------------------------

!  eval_log_primary() exit

return

! **********************************************************************

!  eval_log_primary()

end subroutine eval_log_primary

! **********************************************************************
! **********************************************************************

!  %%% string utilities- editing, parenthesis and quotes

! **********************************************************************
! **********************************************************************

!  replace_substring() edit source lines

subroutine replace_substring( mixed_case_str, lower_case_str, search_str, replace_str, first_idx)

! **********************************************************************

!  replace_substring() interface

! ----------------------------------------------------------------------

!  mixed case string to be printed

character( len= *), intent( inout), optional :: mixed_case_str

!  lower case string to be searched

character( len= *), intent( inout) :: lower_case_str

!  substring to be replaced

character( len= *), intent( in) :: search_str                        

!  string to replace target

character( len= *), intent( in) :: replace_str                       

!  location of first occurance

integer, intent( in) :: first_idx

! **********************************************************************

!  entry: line is a line of Fortran source with (possibly) ?target?

!  exit: line has any ?target? strings replaced

! **********************************************************************

!  replace_substring() local

! ----------------------------------------------------------------------

!  beginning and end of target within lines

   integer :: end_idx

   integer :: search_idx

   integer :: search_len

! **********************************************************************

!  replace_substring() text

! ----------------------------------------------------------------------

continue

! ----------------------------------------------------------------------

!  initialize

   search_idx = first_idx

   search_len = len( search_str)

! ----------------------------------------------------------------------

!  if mixed case is present

   mixed_present: if( present( mixed_case_str) )then

!  replace in both strings

      edit_mixed: do while( search_idx > 0 )

         end_idx = search_idx + search_len

         end_mixed: if( search_idx == 1 )then

            mixed_case_str = replace_str // mixed_case_str( end_idx: )

            lower_case_str = replace_str // lower_case_str( end_idx: )

         elseif( end_idx > len( lower_case_str) )then end_mixed

            mixed_case_str = mixed_case_str( : search_idx - 1) &
                           // replace_str

            lower_case_str = lower_case_str( : search_idx - 1) &
                           // replace_str

         else end_mixed

            mixed_case_str = mixed_case_str( : search_idx - 1) &
                           // replace_str &
                           // mixed_case_str( end_idx: )

            lower_case_str = lower_case_str( : search_idx - 1) &
                           // replace_str &
                           // lower_case_str( end_idx: )

         endif end_mixed

         search_idx = index( lower_case_str, search_str)

      enddo edit_mixed

! ----------------------------------------------------------------------

!  mixed case is not present

   else mixed_present

!  replace in lower case only   

      edit_string: do while( search_idx > 0 )

         end_idx = search_idx + search_len

         end_lower: if( search_idx == 1 )then

            lower_case_str = replace_str // lower_case_str( end_idx: )

         elseif( end_idx > len( lower_case_str) )then end_lower

            lower_case_str = lower_case_str( : search_idx - 1) // replace_str

         else end_lower

            lower_case_str = lower_case_str( : search_idx - 1) &
                           // replace_str &
                           // lower_case_str( end_idx: )

         endif end_lower

         search_idx = index( lower_case_str, search_str)

      enddo edit_string

   endif mixed_present

! ----------------------------------------------------------------------

!  replace_substring() exit

return

! **********************************************************************

!  replace_substring()

end subroutine replace_substring

! **********************************************************************
! **********************************************************************

!  seek_close_paren() true if successfully found matching ()

subroutine seek_close_paren( string, start, match)

! **********************************************************************

!  seek_close_paren() interface

! ----------------------------------------------------------------------

!  the string starting with open parenthesis

character( len= *), intent( in) :: string

!  the index of the open parenthesis

integer, intent( in) :: start

!  the index of the matching close parenthesis

integer, intent( out) :: match

! **********************************************************************

!  seek_close_paren() constants

! ----------------------------------------------------------------------

!  close parenthesis

   character( len= *), parameter :: close_paren = ')'

! **********************************************************************

!  seek_close_paren() local

! ----------------------------------------------------------------------

!  counters and pointers

   integer :: level

   integer :: string_len

! **********************************************************************

!  seek_close_paren() text

continue

! ----------------------------------------------------------------------

!  initialize

   string_len = len_trim( string)

! ----------------------------------------------------------------------

   level = 0

   search: do match = start + 1, string_len

      levels: select case( string( match: match) )

      case( open_paren) levels

         level = level + 1

      case( close_paren) levels

         eureka: if( level == 0 )then

            exit search

         endif eureka

         level = level - 1

      end select levels

   enddo search

! ----------------------------------------------------------------------

!  seek_close_paren() exit

return

! **********************************************************************

!  seek_close_paren()

end subroutine seek_close_paren

! **********************************************************************
! **********************************************************************

!  unquote_string() true if extracts string from between quotes

subroutine unquote_string( quoted_str, unquoted_str, in_len, out_len)

! **********************************************************************

!  unquote_string() interface

! ----------------------------------------------------------------------

!  the quoted string to be unquoted

character( len= *), intent( in) :: quoted_str

!  the unquoted string

character( len= *), intent( out) :: unquoted_str

!  the length of the quoted string

integer, intent( out) :: in_len

!  the length of the unquoted string

integer, intent( out) :: out_len

! **********************************************************************

!  unquote_string() local

! ----------------------------------------------------------------------

!  which quote is to be used

   character( len= 1) :: quote

! **********************************************************************

!  unquote_string() text

continue

! ----------------------------------------------------------------------

!  which quote is the first quote (if either)

   which_quote: select case( quoted_str( 1: 1) )

! ----------------------------------------------------------------------

!  string delimited by single quote

   case( single_quote) which_quote

      quote = single_quote

! ----------------------------------------------------------------------

!  string delimited by double quote

   case( double_quote) which_quote

      quote = double_quote

! ----------------------------------------------------------------------

!  string delimited by neither quote- nothing to do

   case default which_quote

      in_len = 0

      out_len = 0

      unquoted_str = null_string

      return

   end select which_quote

! ----------------------------------------------------------------------

!  initialize scan loop

   in_len = 2
   out_len = 1

   unquoted_str = blank

!  scan thru the quoted string

   scan_string: do while( in_len <= len_trim( quoted_str) )

! ----------------------------------------------------------------------

!  if find one matching quote

      next_char: if( quoted_str( in_len: in_len) == quote )then

!  check next character

         in_len = in_len + 1

!  check for a pair of quotes

         next_quote: if( quoted_str( in_len: in_len) == quote )then

            unquoted_str( out_len: out_len) = quoted_str( in_len: in_len)

            in_len = in_len + 1

            out_len = out_len + 1

         else next_quote

            exit scan_string

         endif next_quote

! ----------------------------------------------------------------------

!  character is not a matching quote

      else next_char

         unquoted_str( out_len: out_len) = quoted_str( in_len: in_len)

         in_len = in_len + 1

         out_len = out_len + 1

      endif next_char

   enddo scan_string

! ----------------------------------------------------------------------

!  unquote_string() exit

return

! **********************************************************************

!  unquote_string()

end subroutine unquote_string

! **********************************************************************
! **********************************************************************

!  seek_match_quote() true if extracts string from between quotes

subroutine seek_match_quote( quoted_str, start_char, in_len)

! **********************************************************************

!  seek_match_quote() interface

! ----------------------------------------------------------------------

!  the quoted string to be unquoted

character( len= *), intent( in) :: quoted_str

!  the length of the quoted string

integer, intent( in) :: start_char

!  the length of the quoted string

integer, intent( out) :: in_len

! **********************************************************************

!  seek_match_quote() local

! ----------------------------------------------------------------------

!  which quote is to be used

   character( len= 1) :: quote

! **********************************************************************

!  seek_match_quote() text

continue

! ----------------------------------------------------------------------

!  which quote is the first quote (if either)

   which_quote: select case( quoted_str( start_char: start_char) )

! ----------------------------------------------------------------------

!  string delimited by single quote

   case( single_quote) which_quote

      quote = single_quote

! ----------------------------------------------------------------------

!  string delimited by double quote

   case( double_quote) which_quote

      quote = double_quote

! ----------------------------------------------------------------------

!  string delimited by neither quote- nothing to do

   case default which_quote

      in_len = 0

      return

   end select which_quote

! ----------------------------------------------------------------------

!  initialize scan loop

   in_len = start_char + 1

!  scan thru the quoted string

   scan_string: do while( in_len <= len_trim( quoted_str) )

! ----------------------------------------------------------------------

!  if find one matching quote

      next_char: if( quoted_str( in_len: in_len) == quote )then

!  check next character

         in_len = in_len + 1

!  check for a pair of quotes

         next_quote: if( quoted_str( in_len: in_len) == quote )then

            in_len = in_len + 1

         else next_quote

            exit scan_string

         endif next_quote

! ----------------------------------------------------------------------

!  character is not a matching quote

      else next_char

         in_len = in_len + 1

      endif next_char

   enddo scan_string

! ----------------------------------------------------------------------

!  seek_match_quote() exit

return

! **********************************************************************

!  seek_match_quote()

end subroutine seek_match_quote

! **********************************************************************
! **********************************************************************

!  coco

! $Id: coco.f90,v 1.30 2007/06/25 19:08:22 dan Exp dan $
! **********************************************************************

end program coco
