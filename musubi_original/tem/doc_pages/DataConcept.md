title: Data Concept

# Sketch of the data concept deployed in treelm

Illustrated with a small example to show the tree numbering, and distribution
in 2D. An overview image of the example is given in 2DIllustration.svg.
In the middle is a obstacle (4 elements on the finest level), and
the domain is distributed to 6 processes.

A breadth first numbering scheme of all nodes in the complete tree is
used, with the positive half of a 64 bit integer, this leaves
us with 20 bits for the finest representable level in each
direction. (1 bit is used by the upper part of the tree).
A second 8 byte integer is used to provide a pointer to possibly
further informations about the element.
The integer in this second entry provides the index of this element
in the list of all elements with special properties, like boundary conditions,
deformations or attached material properties. If it is 0, there
is no additional information.

With this hierarchical scheme of storing the data it is possible to directly
access each information of an element in a distributed manner, but store only
a minimal set of informations for each element.

Neighborhoods are implicitly given by the tree topology, or if necessary has
to specified separately for the elements in the additional property
informations. The logic here is, that the neighborhood needs to be reconstructed
anyway in the case of elements on different levels. And with the help of putting
the neighborhood informations separately if necessary it is possible to support
arbitrary unstructured meshes, with any kind of connectivity between the
elements.

To identify which elements are located on which processes, the
outermost left and right treeIDs of the process local leaves are
gathered on each process (allgather).

In the given example this would look like this:

```fortran
treeID_left =  [  5, 34, 45, 13, 15, 72]
treeID_right = [ 33, 10, 12, 60, 71, 20]
```

From this information, the boundary of each partition can be derived on each
level of the tree:
The parent of a any node is given by the following integer operation:

```fortran
parent = (treeID-1) / 4 ! 2D
parent = (treeID-1) / 8 ! 3D
```

For example the node with treeID 8 has to be found for element 3 (tree ID = 7)
in all the partitions.

* Get the upper boundary for the level, the neighbor is located in:

```fortran
  level = 0
  ub = 1 ! highest tree ID of the searched level
  n = treeID ! 8 in this example
  do while (n > 0)
    level = level + 1
    ub = ub + 4**level
    n = (n - 1) / 4 
  end do
```

* Now get the parents of the bounding leaves on each partition, until
  they are at least on the same level, as the searched tree ID:

```fortran
  do iRank=0,nRanks
    left(iRank) = treeID_left(iRank)
    do while (left(iRank) > ub)
      left(iRank) = (left(iRank) - 1) / 4
    end do

    right(iRank) = treeID_right(iRank)
    do while (right(iRank) > ub)
      right(iRank) = (right(iRank) - 1) / 4
    end do
  end do
  
  !! In the example, this results (for treeID 8) in
  !! left  = [ 5,  8, 11, 13, 15, 17 ]
  !! right = [ 8, 10, 12, 14, 17, 20 ]
  !! That is, it is now known which partition contributes
  !! to which interval on the given level:
  !! rank 0:  5 -  8
  !! rank 1:  8 - 10
  !! rank 2: 11 - 12
  !! rank 3: 13 - 14
  !! rank 4: 15 - 17
  !! rank 5: 17 - 20
  
  !! When looking for the construction of node 8 by its childs it is therefore
  !! necessary to ask rank 0 and rank 1 for their contribution.
```

With this only the local elemental data necessary communication structures need
to be stored on each process, the tree does not have to be actually stored.

The required data on disk is just the [[treelmesh_type:treeID]] and an additional bit field, to
provide the possibility to indicate further properties for each element.
For example boundary conditions would be such a property.
Thus there are just 16 byte per element + additional informations where needed.
To find additional elemental data, where needed, it is necessary to build a
prefix sum, for each used property in order to find the local offset, as direct
links for each of the possible 64 properties would be too large. As an allgather
is required anyway, the much cheaper prefix sum shouldn't be a problem anyway.

A comment on visualization output: it should be straight forward to allow the
output of the mesh down to a specified tree level only, allowing drastically
reduced output file sizes.
