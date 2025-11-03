title: Example for the stencil construction

To build the stencil for a given element, all required
[[treelmesh_type:treeID]]s can be obtained with
the following procedure:
First the integer coordinate tuple \((x,y,z,L)\) of the
[[treelmesh_type:treeID]]s, for which the stencil
has to be obtained, is computed.
Then for each stencil element the
[[treelmesh_type:treeID]] is obtained by the
following steps:

- Add the stencil offset to the coordinate of the central element to obtain
  the coordinate \((x + s_x, y + s_y, z + s_z, L)\) for the stencil element

- Convert the obtained coordinate back to a
  [[treelmesh_type:treeID]].

To illustrate this procedure with a specific example, lets consider a two
dimensional mesh.
We would like to find the right neighbor of the element with
[[treelmesh_type:treeID]] \( 40 \).
An illustration of the geometrical layout for this example is given in the
following Figure.
As we are looking at a two dimensional mesh, there are \( 4^L \) elements on
each refinement level \( L \).

![Illustration of the right neighbor identification](|media|/find_rightNeighbor.png)

The first information, that is needed, is the level of the given
[[treelmesh_type:treeID]].
To find this, and the Morton index on that level, we subtract the elements
on coarse levels, until the reduced index is smaller than the next number of
elements to subtract:

- Level 0: \(40 - 4^0 = 39\)
- Level 1: \(39 - 4^1 = 35\)
- Level 2: \(35 - 4^2 = 19\)
- Level 3: \(4^3 > 19\)
- [[treelmesh_type:treeID]] \(40\) is located
  on Level 3 and has a Morton index of \(19\).

The binary representation of \(19\) is \(010011_b\).
Due to the bit-interleaving rule, odd bits correspond to the X coordinate and
even bits correspond to the Y coordinate.
Thus we obtain for X \(101_b = 5\) and for Y \(001_b = 1\).
Now we found the coordinate tuple \((x, y, L) = (5, 1, 3)\) of
[[treelmesh_type:treeID] \(40\) and can
simply add the offset \((1, 0, 0)\), that represents the right neighbor.
By this we obtain the coordinate of this neighbor to be \((6, 1, 3)\),
resulting in a binary representation of \((110_b, 001_b)\).
After interleaving the bits, the Morton index of the right neighbor is found
to be \(010110_b = 22\)
Finally, by adding the offset, we obtain the
[[treelmesh_type:treeID]] of the right neighbor
as \(22 + 4^0 + 4^1 + 4^2 = 43\).
Note, that there is no restriction on the offsets, that are to be applied in
the neighborhood search, and due to the periodic universe all offset
locations are well defined.
