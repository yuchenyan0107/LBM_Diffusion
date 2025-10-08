title: Dynamic data structures

During the processing of the mesh to create the solver data structures,
common tasks arise.
These include among others
- a unique collection of elements in lists, i.e. without duplicates,
- a localization of a stored element,
- sorting of a list of records.
As Fortran does not intrinsically include dynamically growing data
structures, they first have to be introduced manually.

Therefore the modules [[tem_dyn_array_module]] and [[tem_grow_array_module]]
have been introduced. These modules are based on generic CoCo macros stored in
arrayMacros.inc.

The definition is depicted in listing lst:tem_dyn_array where macros
are used for the label of a dynamic type `?tname?`and the used data type
`?tstring?`.
@todo What shall be listed in lst:tem_dyn_array ??

The introduction of macros allow for a generic definition for all data types of
the container, whereas the compiled source code contains the definitions for
each required data type.

# Dynamic Array
The dynamic array is a data structure to store sets of unique elements in a
sorted manner. Elements can be added dynamically, and the dynamic array takes
care of allocating the space needed to store the elements as well as checking
them for uniqueness. Also sorting is done without user interaction.

## Initialization
The initialization is not needed in all cases. Usually, the dynamic array can be
used out of the box. In case you know the total number of elements in advance,
you can provide this information via the init routine. With this hint, the data
structure is able to allocate enough space for the expected number of elements
and thus saves some resizing steps.
When `init` is called, it simply allocates the container with the given options
and resets the number of records.

    type(dyn_intArray_type) :: dynIntArray
    call init( me = dynIntArray, length = nTotalElements )

## Appending elements
To append an element, simply call `append`. As a dynamic array only contains
unique elements, the routine checks for the existence of the entry in the
list and, subsequently, the sufficiency of the size of the container to
incorporate another entry. In case the size does not suffice, the array is
expanded. Information about if (`wasAdded`) and where (`pos`) the entry was
stored is returned to the caller.

    integer :: newPos
    logical :: wasAdded
    call append( me = dynIntArray, val = 1, length = 5, pos = newPos, wasAdded = wasAdded )

Also appending more than one value is possible using the vectorial append
routine:

    integer :: newPos(:)
    logical :: wasAdded(:)
    call append( me = dynIntArray, val = ( 2, 3, 4), length = 5, pos = newPos, wasAdded = wasAdded )

## Expanding a dynamic array
An uninitialized dynamic array starts with the size of 0, but also an
initialized array can become to small over time. If elements should be added,
but the dynamic array doesn't have any free spots left, `expand` is called to
increase the array's size.

If `expand` was called and a new length was provided via the `length` argument,
this length is taken to increase the array. If no length was provided, the array
is simply doubled, but at least increased by [[env_module:minLength]].

Another thing influencing the size is the amount of elements that are added. The
array is always increased to afterwards contain all new elements.

But the array is only increased up to the maximum value of the index,
which means that an dynamic array can only contain up to `huge(me%nVals)`
elements, thus \(2^31\) elements.


## Retrieving the position of an element
You might want to retrieve the position of an element within a list. In order to
achieve this efficiently, the list with \(n\) entries must be ordered so we can
rely on the binary search, which performs with \(\mathcal{O}( \log n )\).

Hence we introduce a derived data type which allows an efficient and simple
element adding and retrieving process, as well as a sorted order of the
entries.
[[tem_dyn_array_module]] provides the type definitions and all the related
procedures to act on this data type.
The data type includes the size of the container within, the amount of its
actual entries, a flag for uniqueness and the container array itself as well as
a conforming array for the sorted retrieval of the entries in the value array.

## Truncating a dynamic array
See [Truncating a growing array](#Truncating_a_growing_array).

## Emptying a dynamic array
See [Emptying a growing array](#Emptying_a_growing_array).

## Destroying a Dynamic array
See [Destroying a growing array](#Destroying_a_growing_array).

# Growing Array
A growing array is a data structure to store arbitrary sets of elements. The
growing array does not check for uniqueness nor does it sort it's elements. Thus
it is similar to the common fortran array except that it can grow dynamically
(at least that's how it looks like).

## Initialization
The initialization of the growing array the same applies as it does for the
dynamic array; an initialization is only needed, when th ethe total number
of elements in known in advance.

    type(grw_intArray_type) :: grwIntArray
    call init( me = grwIntArray, length = nTotalElements )

## Appending elements
Elements can be appended in two different ways. The first one adds a single
element to the grwoing array:

    call append( me = grwIntArray, val = 1 )

The second routine takes an array of values:

    call append( me = grwIntArray, val = ( 2, 3, 4 ) )

After these two calls, the grwoing array contains `( 1, 2, 3, 4 )`.

## Placing elements at a given position
When elements should end up at a speciic position, `append` is not the routine
of choice, as it appends the new elements to the end of the array. Use `placeAt`
instead. `placeAt` needs to be called with the additional argument `pos`, which
states the position the element should be placed at.

Regarding the example from "Appending elements", the current array looks like
`( 1, 2, 3, 4 )`. After

    call placeAt( me = grwIntArray, val = 5, pos = 3 )

the array wil look like `( 1, 2, 5, 4 )`. There is also a vectorial version
available:

    call placeAt( me = grwIntArray, val = ( 5, 6, 7 ), pos = 3 )

which results in `( 1, 2, 5, 6, 7 )`.

## Expanding the growing array
The array starts with a size of 0. When a new element is appended, the array
checks whether the new element fits into it. If not, the array is expanded.
If the caller doesn't explicitly provide a length, the array will be doubled
(whereas the minimum expansion length is [[env_module:minLength]]).
Regarding the example from the paragraph "Appending elements", the array
contains 4 elements and the array is also 4 elements long. When a fifth element
is added, the array will be resized, so that `size(grwIntArray%Val)` will return
8, even when `grwIntArray%nVals` is 5.

In cases where the final size of the array is known or the size of the batch of
elements that will be added, this size can be passed along with calls to append
or placeAt using the length argument. However, passing the expand length doesn't
necessarily mean that it is used. If it is smaller than the position requested
by placeAt is even outside the expanded array or the expanded array would exceed
it's maximum size, other expand lengths will be used.

Liek the dynamic array, the maximum size for the growing array is also
`huge(me%nVals)`, which denotes to \(2^31\) elements.

## Truncating the growing array
With the rather involved logic to expand a growing array, it can happen that
the final array has much more element slots available that real elements (which
means: `me%nVals` << `me%ContainerSize`). To free the surplus, the `truncate`
routine shrinks the array to the actually needed size.

Regarding the array `( 1, 2, 5, 6, 7, x, x, x )` with x indicating the available
slots. After

    call truncate( me = grwIntArray )

the array will be `( 1, 2, 5, 6, 7 )` and therefore
`me%nVals` == `me%ContainerSize`.

## Emtying a growing array
If a growing array should be recycled, but it's elements are not needed anymore,
the array can be emptied by calling

    call empty( me = grwIntArray )

This actually doesn't change the array's container. It will stay untouched,
which means: If it was 8 slots long before calling `empty`, it is also 8 slots
long after the call. The growing array just forgot that the container contains
anything useful.

## Destroying a growing array
In case a growing array is not needed anymore, the memory consumed by the
elements with the array and the growing array itself can be freed by calling

    call destroy( me = grwIntArray )

# Two-dimensiona Growing Array
There is also a two-dimensional growing array, which actually is an array of
growing arrays. The width of the two-dimensional growing array can only be set
during initialization of the array, whereas the sizes of each array can be set
individually and also after the initialization is already done (or not, as a
growing array doesnot necessarily has to be initialized before use).

The following paragraphs do only state the differences to the common,
one-dimensional array.

## Initialization
Initialization is similar to the common growing array, with the difference that
the width for the first dimension is mandatory (length stays optional):

    call init( me = growing2dIntArray, width = 3, length = 3 )

## Appending elements
To append an element, the index of the first dimension has to be provided:

    call append( me = grw2dIntArray, val = 1, pos1 = 2 )

After initializing an array with 3x3 elements and calling the above example, the
array looks like:

    ( ( _, _, _ ), ( 1, _, _ ), ( _, _, _ ) )

Providing length and retrieving the position where the element was added is also
possible for the two-dimensional version:

    integer :: addedPos
    call append( me = grw2dIntArray, val = 1, pos1 = 2, length = 3, pos2 = addedPos )

Also vectorial access is possible. But in contrast to the one-dimensional
version, not all elements are added to one growing array, but to the first free
index for all growing arrays. Consider the array from the exmaple above:

    ( ( _, _, _ ), ( 1, _, _ ), ( _, _, _ ) )

Calling the vectorial append

    call append( me = grw2dIntArray, val = ( 4, 5, 6 ) )

will result in

    ( ( _, 4, _ ), ( 1, 5, _ ), ( _, 6, _ ) )

because the element `1` in the second array blocks the first index, thus the
second index is the first one that is free for all arrays.

## Placing elements at a given position
When placing elements at a requested position, the scalar approach behaves as
expected. The vectorial access is, like for appending values, a little
different. Taking the array from the example above:

    ( ( _, 4, _ ), ( 1, 5, _ ), ( _, 6, _ ) )

calling

    call placeAt( me = grw2dIntArray, val = ( 7, 8, 9 ), pos2 = 2 )

will result in

    ( ( _, 7, _ ), ( 1, 8, _ ), ( _, 9, _ ) )
