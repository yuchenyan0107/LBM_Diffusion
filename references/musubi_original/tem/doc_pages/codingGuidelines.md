title: Coding Guidelines


## Fortran

Writing code is a complex task. But it is not only complex in the case of the
problem it tries to solve, but also in the case of the code itself.
Writing good code is more than just writing functional code.
Code quality is related to the correctness of the results as well as to the
maintainability of the code.
Keep in mind, that you are not the only one working with the code.
And even if you are currently the only one working with your code this does not
imply that this will continue.
So try to be a nice guy and help others to understand your work.
The rule *It was hard to write, so it has to be hard to read* is not the best
principle when it comes to code design.
Also remember, that you might need to read the code yourself after some time,
and then it usually is helpful to have a code as readable as possible.


### Files

Each source file should contain only one top level program unit, and should be
named after the program unit it contains.
Mostly this means one module per file.
As we are using the free source form, the file extension is *f90*.
When a preprecessor like [CoCo](http://www.daniellnagle.com/coco.html) is used,
*fpp* has to be used as file extension.
This is a technical necessity, as waf will process the files based on these
suffixes.


### Coco usage

If you need to reimplement functionality without coco, you may use it
(except for in Aotus). Otherwise please avoid its usage. Specifically, do not
use it to short-cut coding only because you have similar looking code parts
that need to be replicated multiple times. Thus, growing and dynamic arrays
for example are ok, as you would need to implement the same stuff for
different datatypes. Similarly, array access via a macro are OK,
as you would need to implement basically all code parts using the array
with the different access patterns. NOT OK is replacing interfaces of
routines, which you anyway have to implement a new. Those interfaces are
just partial code fractions and replacing them with coco templates is simply
shortcutting the code.


### Comments

Comments should describe the purpose of the code and what a specific section
is meant to achieve.
Comments describing what a code does are usually of little use, as this already
presented by the code itself.

We use [FORD](https://github.com/cmacmackin/ford) to automatically generate code
documentation from the comments in the code.
Please use the markers from FORD to properly fill this automatic
documentation.
Generally a FORD comment is started by the `!!`.
If you want to put a comment in front of the entity to add the documentation to,
you have to use `!>` to start the FORD comment.

Even though FORD documentation usually starts with '!!‘, there are several other
characters which have a specific meaning in the second place (like ‚!*‘),
you can even define them yourself in the FORD configuration.
So for those formatting comments it is better to always leave a space
in the second place. Also, please indent the comments along with the code.

For inline code you need to start the a block with three ticks, which need to be
immediately attached to the exclamation marks at the beginning of the line.
You may also specify a language for the syntax highlighting like this:
```fortran
!!```lua
!!  somelua = {}
!!```
```

It is also possible to include other files within code blocks by means of the
[markdown include plugin](https://github.com/cmacmackin/markdown-include).

In our code, we use `! ************* !` to separate subroutines and
functions and `! ———— !` to separate input variables and local variables. These
separators are also indented along with the code it is separating. See also the
[structure examples](#structure).


### Control structures

If you have deeply nested control structures (Do-loops, ifs, cases), use block
names to differentiate them.
This also allows you to directly address outer loops from inner loops, and
thus enables you to exit a complete nested loop block early.
These names allow the compiler to check the constructs on consistency and at the
same times allows the reader an easier identification of the individual part
associations.


### Module structure

At the very beginning of each file containing a module, there has to be a
comment describing the module, its purpose, as well as its usage and perhaps its
pitfalls.

Other modules are imported using the `use` keyword.
To have the chance to identify the source of a foreign subroutine, function, or
variable, modules should not be used as a whole, but only the parts that are
really needed in the current module.
This is achieved by using the `only` keyword, followed by the needed
identifiers.

Each module has to contain an `implicit none` statement to ensure that there is
no implicit type definition.
Any code that is not covered by an `implicit none` might be considered as
erroneous as implicit typing might lead to unexpected results.

Additionally, each module should be made private by default using a `private`
statement to  control which parts of the code are accessible by other modules.
Especially this restricts the accessible entities of the module, such that
those entities used from other modules are not available to using modules.
To export symbols, explicitly state them to be public.

In summary, the module dependency is kept as explicit as possible, and,
therefore, the code should be least surprising to the reader.
This also helps the compiler to optimize the code.

This leads to the following basic structure for a module:

```fortran
!> Here is the module's short description.
!!
!! And here comes the module's longer description, which can also span over
!! multiple lines.
!! If it gets too long it might be better to separate the details into a
!! separate page and link that page here.
module example_module

  use anotherModule, only: some,           &
    &                      thing,          &
    &                      somethingOther, &
    &                      things

  implicit none

  private

  integer, parameter :: anythingToBePublic = 1

  public :: anythingToBePublic, action


contains


  function action()
    call somethingOther()
  end function action

end module example_module
```


### Subroutines and functions

Always declare an argument with an `intent` attribute.
The only exception is when the argument is a pointer or a procedure.
A function should have only `intent(in)` arguments.
Use a subroutine if more than one argument needs to have an `intent(out)` or any
argument needs to have an `intent(inout)` attribute.
Also, beware of functions with respect to OpenMP parallelism, as there might
arise problems here, if arrays are to be assigned with the result of the
function.

For better readability, always use named arguments to call a function or
subroutine, unless there is only one mandatory argument. This rule maintains a
clearer image on the arguments that are common across all calls. In cases where
the meaning of the single argument provided by a function or subroutine
call is not immediately obvious, also use the named argument calling convention.

Instead of `call tem_time_load(me%min, conf, 'min', thandle)`, the call should
look as follows:

```fortran
call tem_time_load( me     = me%min, &
  &                 conf   = conf,   &
  &                 key    = 'min',  &
  &                 parent = thandle )
```

When you are passing an array argument to a procedure, the dummy argument
usually should be an assumed shape array. Only if this is not possible, other
dummy array argument specifications may be used.
Memory management of assumed shape arrays is automatically handled by the
procedure.
By using the assumed shape definition, you allow the compiler to make use of
array descriptors, which avoid copy operations, even if the calling side passes
in a fragmented array.
You also inherit a matching array description within the routine, which enables
bounds checks by the compiler.


#### <a name="structure"></a>Structure

Every subroutine and every function has a preceeding and a trailing comment line
consisting of stars from the routine's indentation up to column 80 minus the
indentation with a closing exclamation mark:

```fortran
  ! ************************************************************************ !
  subroutine do_something(some, arguments)
    integer, intent(in) :: some, arguments
    integer :: i, nDiagonals, diag_off

    some = code_here(arguments)

  end subroutine do_something
  ! ************************************************************************ !
```

Separate the procedures by two empty lines.


#### Argument and local variable list

The arguments, as well as the local variables, are packed into lines of dashes,
also starting (and ending) at the current indentation:

```fortran
  ! ************************************************************************ !
  subroutine do_something(some, arguments)
    ! -------------------------------------------------------------------- !
    integer, intent(in) :: some, arguments
    ! -------------------------------------------------------------------- !
    integer :: i, nDiagonals, diag_off
    ! -------------------------------------------------------------------- !

    some = code_here(arguments)

  end subroutine do_something
  ! ************************************************************************ !
```

For each argument there is a describing (FORD) comment.
Variable definitions should also be accompanied by some describing comment, if
the name itself is not explaining itself already.
FORD will not make use of local variable definition comments, so do not use
the FORD syntax for these.
Variables with their comments should be separated by blank lines from each
other to have a better indication where the definition of a certain variable is.


```fortran
  ! ************************************************************************ !
  subroutine do_something(some, arguments)
    ! -------------------------------------------------------------------- !
    !> Comments for arguments are mandatory. In most cases arguments have some
    !! constraints, or they are changes during method execution or... So be a
    !! nice guy and provide the user of your code with the information he needs
    !! to get his and your code running without making mistakes.
    integer, intent(in) :: some

    !> Making mistakes is a good keyword. Provide information about pitfalls
    !! one has to take in mind when using your method or passing an argument.
    integer, intent(in) :: arguments
    ! -------------------------------------------------------------------- !
    integer :: i ! (obvious loop index, not need for commenting)

    ! The number of diagonals. This is a borderline case as the variable name
    ! is quite explaining. One can argue that n is not to be read as number of
    ! by all means, though we typically use it in this way.
    integer :: nDiagonals

    ! Here again we have a partly self describing variable name. But in this
    ! case it is just more obvious. One who knows the code does also know that
    ! off stands for offset. But comments are not only for those who already
    ! know what the code does.
    integer :: diag_off
    ! -------------------------------------------------------------------- !

    some = code_here(arguments)

  end subroutine do_something
  ! ************************************************************************ !
```


### Formatting

Use the Fortran 90 free form syntax.
Use spaces and blank lines where appropriate to format your code to improve
readability.


#### Lines

The line width has to be confined to 80 characters.
Short and simple lines are easier to parse and understand than longer ones.
The only valid exception here might be overly long CoCo expressions, which need
to be on a single line.
Use two empty lines before and after contains statements.


##### Continuation

If a line exceeds the maximum line length, it has to be continued using the
Fortran line continuation syntax.
Thus, each continued line ends with the ampersand character.
We also begin the proceeding line an indented ampersand.
The leading, as well as the trailing ampersands, have to be aligned.
The leading ampersands are also indented using the common indentation width of
two spaces.
Closing brackets have to be aligned with the trailing ampersands.
Thereby, the code becomes an easily recognizable visible block.

```fortran
call a_method( with     = a,       &
  &            very     = long,    &
  &            argument = list,    &
  &            that     = exceeds, &
  &            the      = maximum, &
  &            line     = length   )
```

For line continuation, use operators: `//`, `%`, `+`, `-`, `*`, `/`
etc at the beginning of the continued line.
If you need to break a long formula, start the continued lines with operators
to make this continuation more obvious on the next line.

Also, if you need to break an array segment specification, start the continued
line with the `:` to make this immediately visible.

When you need to break a deep derived data-type addressing put the `%` below 
the last % of the preceeding line. If it still doesn't fit, split the preceeding
line again. In case the preceeding line doesn't contain a `%`anymore, indent it
by two starting from the variable name.
Although, you should in addition think about, why this deep access is needed
at all, as the details of derived data types usually should not be of a
concern to the outer procedures.

```fortran
call atl_preprocess_modg_kernel(                          &
  &    equation           = equation,                     &
  &    statedata          = statedata_list(currentLevel), &
  &    mesh               = mesh_list(currentLevel),      &
  &    boundary           = boundary_list(currentLevel),  &
  &    scheme             = scheme_list(currentLevel),    &
  &    material           = material_list(currentLevel),  &
  &    poly_proj_material = poly_proj_list(               &
  &                           material_list(currentLevel) &
  &                             %poly_proj_pos)           )

call atl_modg_modalVolToModalFace(                         &
  &    nElems_fluid = mesh_list(currentLevel)              &
  &                     %descriptor                        &
  &                     %elem                              &
  &                     %nElems(eT_fluid),                 &
  &    length       = mesh_list(currentLevel)%length,      &
  &    volState     = statedata_list(currentLevel)%state,  &
  &    faceRep      = facedata_list(currentLevel)%faceRep, &
  &    nScalars     = equation%varSys%someMoreLevel        &
  &                                  %nScalars,            &
  &    orlikethis   = orlikethis%varSys                    &
  &                             %someMoreLevel             &
  &                             %nScalars,                 &
  &    nDerivatives = equation%nDerivatives,               &
  &    modg         = scheme_list(currentLevel)%modg       )
```

Similarly, to break a string put the `//` at the begining of the continued line.
Do not split lexical tokens across line continuations!
Do not split a character string across line continuations, as this does not
fit to our indentation rules.

Also, please put line continuation into alignment.
For one continued line the & should be placed into the same column.

Please note, that there is a limit to the continued lines in the Fortran
standard (Fortran 90: up to 39 continuation lines).
Even so Fortran 2003 allows more continuation lines, we should try to stick to
the limit of 39, beyond that it also gets hard to read the code.


#### Indentation

Indent blocks by two spaces.
Tabs are not allowed for indentation, and actually not at all in the complete
source file.
Be aware, that this is technically enforced by mercurial, and you will not be
allowed to push Fortran source files, that contain any tabs.
Comments should be indented with the code.

To improve readability and support visual block structures, indentation can be
aligned with the preceeding construct. For example, if you linebreak a nested
array access, indent the `%` two spaces deeper than the preceeding variable
name (see the example above with `nElems_fluid`).

Those visual blocks can also be created with the named argument lists when the
call is split onto several lines. In the above example, `nElems_fluid` is two
spaces _deeper_ than the subroutine name.

#### <a name="spaces"></a>Spaces

* Keywords and operators are wrapped into spaces.
  The statement `if((a+b)>=(c+d))then`, therefore, should be written as
  `if ((a + b) >= (c + d)) then`.
* For nested brackets it also might be useful to add spaces, to more clearly
  indicate the grouping.
  Always use wider spacing for the brackets on outer levels.
  The example above would become `if ( (a + b) >= (c + d) ) then`, though in
  this case the expression is that simple, that the additional spacing is not
  really necessary.
* End statements are separated from the block type specifier by one space,
  e.g. `end if` instead of `endif`.
  The latter notation is a backward compability and not available for newer
  constructs like select or where. The notation with a separating space is
  therefore more consistent.
* Arguments in an argument list are at least separated by comma and space:
  `call do_something(a, b, even_c)`, not `call do_something(a,b,even_c)`.
  Usually you should use keywords and put each argument on a separate line.
* For select case, do not indent each case but leave a space after each case
  keyword. Example:
```fortran
select case (dim)
case (1)
  ! do something if dim = 1
case (2)
  ! do something if dim = 2
case default
  ! do something if dim is not 1 or 2
end select

select case (scheme_list(currentLevel)%scheme)
case (atl_modg_scheme_prp) ! MODG kernel
  select case ( trim(equation%eq_kind) )
  case ( 'maxwell', 'maxwelldivcorrection', 'euler', 'navier_stokes', &
    &    'filtered_navier_stokes', 'acoustic', 'heat', 'lineareuler', &
    &    'loclineuler'                                                )
  end select
end select
```

### Naming

Give all entities you create a meaningful name.
The best case is that someone does not need to read an explanatory comment to
understand what a routine does, a function returns or a variable contains
because the name already provides this information.

Do not define variables with single character, not even for loop counters.
Use either `ii`, `jj`, `kk` or, for better readability, might be some indication
what the counter is actually iterating over.
For example `iX`, `iY` and `iZ` could provide more explanation to the reader.
This is especially true with many nested or different loops.
For other variables there is usually a descriptive name that can be used.
For example `blockSize` is just more expressive than `s`.
Use `i` and `n` as prefixes, where `i*` indicates a counter,
while `n` indicates a count, like the size of an array.

Public entities in APES modules have to obey some more naming rules, which allow
us to avoid name clashes and easy identification of the entities.
We use a three letter prefix for each project, e.g. `tem` for treelm.
The next part of the name is the module or feature name, e.g. `tracking`.
Finally, the actual description of the entity is appended.
Thus, a valid and meaningful name for a routine is, e.g.,
`tem_tracking_dumpAsciiTransient`.
Please notice the [camel casing](https://en.wikipedia.org/wiki/CamelCase) we
used for the description part of the name which increases readability of
identifiers.
Additionally, for subroutines you should use some verb, as above dumping
describes, what the routine does.
Functions on the other hand should describe the result they provide, for
example `tem_parentof` or `tem_coordOfID`.


#### Magic numbers

Avoid magic numbers! Use constants instead.
`allocate( meshInfoList( 6 ))` is not only less expressive than

```fortran
!> number of mesh variables
integer, parameter :: nMeshVars = 6

allocate( meshInfoList( nMeshVars ))
```

but also easier to change, if you ever happen to need to modify `nMeshVars`.
With literal constants, you need to replace all matching occurences, which might
be difficult if you need the same number with different meanings in different
parts of the code.
While with the named constant, there is just a single place, where you need to
change the code.
Use Magic numbers only if it is local or self-explaining.

### Types

Do not declare variables of type real without specifying a kind:
`real(kind=rk) :: realValue`.
The common real kind is available from the env_module and has to be imported
using `use env_module, only: rk`.
For special cases, there are more kinds available in env_module, so there should
not be a need to define a real kind on your own.
Please note, that not all kinds might be provided by the compiler you use.

### Error handling

Use Errcode from aotus' aot_get_val routines to handle loading errors.
If you don’t pass in ErrCode to the aot_get_val, the library will abort
and each process will write an error
message, which is probably not what you want.
If you pass ErrCode, you are responsible for treating errors yourself,
so please take care of it. You can get
the Lua error message in ErrString.

### Operators

Use the new syntax for operators:

|| New || Old ||
|| --- || --- ||
|| == ||.EQ. ||
|| /= ||.NE. ||
|| \> ||.GT. ||
|| <  ||.LT. ||
|| >= ||.GE. ||
|| <= ||.LE. ||


## Python

See [PEP 8 -- Style Guide for Python Code](http://legacy.python.org/dev/peps/pep-0008/).

For docstrings, use the [numpy format](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt)


## Lua

Lua has at's [own guidelines](http://lua-users.org/wiki/LuaStyleGuide), which
are very similar to Python's.

### Indentation

Indent table contents by two spaces. Put the opening bracket on the line with
the table name and the closing bracket at the same indentation level as the
table name. When the content inside the brackets is rather short, e.g. a vector
or one named component, the content as well as the opening and closing brackets
can stay on one line. If there is more than one named component, split it into
several lines.

``` lua
table = {
  comp = { 1.0, 2.0, 3.0 },
  another_table = {
    with = 'string',
    some = 4711,
    variables = 42.23
  }
}
```

### Spacing

Must:

* Space after comma.

Should:

* Add a space after the opening and in front of the closing brackets.
* [Like in fortran](#spaces), operators should be separated by spaces, as long
  as they don't disturb readability of formulas (e.g. group multiplications,
  split additions).


### Casing

Lua is, in contrast to Fortran, case sensitive. In the context of APES, Lua is
used to configure the different components of the APES software suite written
in Fortran, which leads to possible breaking points due to different casing
paradigms. To prevent misunderstandings between Lua and Fortran, Fortran code
has to normalize all
[values](https://en.wikipedia.org/wiki/Value_%28computer_science%29) (usually
obtained via [[aot_table_get_val]]) and
[identifiers](https://en.wikipedia.org/wiki/Identifier#In_computer_science)
(the name of variables, tables etc.) from the "Lua domain". To avoid
confusion between Lua statements and normalized Fortran statements, it is
recommended to use `lower_case` notation throughout all code files.
