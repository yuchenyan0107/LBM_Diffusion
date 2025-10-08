title: Treelm File Format

As already introduced in the previous sections, the
[octree mesh](octree.html) is represented by the sparse
storage of the leaf nodes in the computational domain.
The storage of the mesh on disk is basically the sequence of [[treelmesh_type:treeID]]s
building the computational domain.
However, this is in general insufficient to fully describe the mesh and an
additional elemental information is introduced to attach further data to the
elements.
TreElM uses 64-bit signed integers to encode the [[treelmesh_type:treeID]]s for
greatest portability.
In addition, a second 64-bit signed integer is used to build a bit-mask for
properties like boundary conditions, that might be attached to each element.
This data is written in native binary format to disk, resulting in a file
with 16 bytes per element.
Due to this simple format, the data might be accessed on an elemental basis
with Fortran's direct IO, it can be easily converted between big and little
endian representations and allows fast reading and writing.

This binary data is accompanied by a header file with descriptive information
about the mesh, like the total number of elements, the physical extent and
origin of the universal cube and descriptions for the attached properties.
The data in this header is provided in the form of a [Lua](http://www.lua.org/)
script, to allow for a flexible handling of additional data and future evolution
of the format.
With the help of the Aotus library ([[aotus_module]]) data can be easily
retrieved from and written to such Lua scripts.
