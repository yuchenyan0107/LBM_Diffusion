title: Face Details

Details and examples for the Face Module
========================================

This page describes some face configurations as they may appear in the
[[tem_face_module]] in more detail.

Serial face description (no level jumps)
----------------------------------------

In case of a serial execution the only face situation that can occur is:

```
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /      rank 0      /      rank 0      /
   /                  /                  /
  /                  /                  /
  --------------------------------------
```

All fluid-fluid faces are compute faces.


Parallel face descritption (no level jumps)
-------------------------------------------

In the case of parallel execution without local refinements,
the only additional face situation that can occur is:

```
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /      rank 0      /      rank 1      /
   /                  /                  /
  /                  /                  /
  --------------------------------------
```

In this situation the rank that owns the left element is the one that
computes the face.

So on rank 0 we have the following situation:

```
                      compute face
                           /
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /                  /                  /
   /                  /                  /
  /                  /                  /
  --------------------------------------
```

Therefore, rank 1 will send face information in the computation and the
siuation is:

```
                         send face
                           /
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /                  /                  /
   /                  /                  /
  /                  /                  /
  --------------------------------------
```


Parallel face description with level jumps
------------------------------------------

In parallel a set of level-jumb-partition-bnd face combinations may occur.
Below, you can find a collection of thse cases (with one level difference)
and how they are categorized.

### Mesh setup ###

The faces are categorized by dimension-by-dimension level descriptors.
Therefore, it is sufficient to consider the different face combinations in
a single spatial direction. Here, we restrict ourselves without loss of
generalty to the x direction. The situation we want to consider is the
following:

```
       ---------------------------------------
      /        /         /                  /
     /        /         /                  /
    /------------------/                  /
   /        /         /                  /
  /        /         /                  /
  --------------------------------------
```

In the following the numbers in the elements will represent the rank the
element is located on.

### Example Rank Setup 1 ###

```
       ---------------------------------------
      /   1    /    1    /                  /
     /        /         /                  /
    /------------------/        0         /
   /   1    /    1    /                  /
  /        /         /                  /
  --------------------------------------
```

#### Layout on rank 0: ####

```
                       from finer face
                           /
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /                  /                  /           coarse level
   /                  /                  /
  /                  /                  /
  --------------------------------------
       ---------------------------------------
      /        /         /        /         /
     /        /         /        /         /
    /------------------/------------------/           fine level
   /        /         /        /         /
  /        /         /        /         /
  --------------------------------------
                   /
              send face
```

#### Layout on rank 1: ####

```
           neither fromFinerFace nor receive face
                           /
       -------------------- - - - - - - - - -
      /                  /                  /
     /                  /                  /
    /                  /                              coarse level
   /                  /
  /                  /
  -------------------
       ---------------------------------------
      /        /         /        /         /
     /        /         /        /         /
    /------------------/------------------/           fine level
   /        /         /        /         /
  /        /         /        /         /
  --------------------------------------
                   /
           compute + recv faces
```

@note Please notice, that this situation is covered (compute + recv) on
the finest face level of rank 1 although the fine elements are located on a
single remote rank from the perspective of rank 0. This is in contrast to the
behavior of the level descriptor where the information is covered on the
coarsest level as possible.


### Example Rank Setup 2 ###

```
       ---------------------------------------
      /   1    /    1    /                  /
     /        /         /                  /
    /------------------/        0         /
   /   2    /    1    /                  /
  /        /         /                  /
  --------------------------------------
```

#### Layout on rank 0: ####

```
                       fromFinerFace
                           /
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /                  /                  /           coarse level
   /                  /                  /
  /                  /                  /
  --------------------------------------
       ---------------------------------------
      /        /         /        /         /
     /        /         /        /         /
    /------------------/------------------/           fine level
   /        /         /        /         /
  /        /         /        /         /
  --------------------------------------
                   /
                sendFaces
```

#### Layout on rank 1: ####

```
       ---------------------------------------
      /        /         /        /         /
     /        /         /        /         /
    /------------------/------------------/           fine level
   /        /         /        /         /
  /        /         /        /         /
          ------------------------------
                   /
       receive and compute faces
```

#### Layout on rank 2: ####

```
       -                ----                --
      /                  /                  /
     /                  /                  /
    /---------         /                  /           fine level
   /        /         /                  /
  /        /         /                  /
  -------------------------------------
```


### Example Rank Setup 3 ###

```
       ---------------------------------------
      /   1    /    1    /                  /
     /        /         /                  /
    /------------------/        0         /
   /   1    /    2    /                  /
  /        /         /                  /
  --------------------------------------
```

#### Layout on rank 0: ####

```
                      from finer face
                          /
       ---------------------------------------
      /                  /                  /
     /                  /                  /
    /                  /                  /           coarse level
   /                  /                  /
  /                  /                  /
  --------------------------------------
       ---------------------------------------
      /        /         /        /         /
     /        /         /        /         /
    /------------------/------------------/           fine level
   /        /         /        /         /
  /        /         /        /         /
  --------------------------------------
                   /
              send face
```

#### Layout on rank 1: ####

```
       ---------------------------------------
      /        /         /        /         /
     /        /         /        /         /
    /------------------/------------------/           fine level
   /                  /                  /
  /                  /                  /
                  -----
                   /
       receive and compute faces
```

#### Layout on rank 2: ####

```
                        ----
      /                  /                  /
     /                  /                  /
    /------------------/------------------/           fine level
   /        /         /        /         /
  /        /         /        /         /
  --------------------------------------
                   /
       receive and compute faces
```

Parallel face description with multi-level-jumps
------------------------------------------------

Mesh setup (with multiple level difference)

As a simple example we will consider face description for the following mesh
setup.

```

          ---------------------------------------------------
         /     /      /     /     /                        /
        /-----/------/-----/-----/                        /
       /     /      /     /     /                        /
      /-----/------/-----/-----/                        /
     /     /      /     /     /                        /
    /-----/------/-----/-----/                        /
   /     /      /     /     /                        /
  --------------------------------------------------

```

i.e. in our analysis we have to consider the following layers:

```

          ---------------------------------------------------
         /                        /                        /
        /                        /                        /
       /                        /                        /      coarse level
      /                        /                        /
     /                        /                        /
    /                        /                        /
   /                        /                        /
  --------------------------------------------------
          ---------------------------------------------------
         /            /           /           /            /
        /            /           /           /            /
       /            /           /           /            /
      /------------/-----------/-----------/------------/   intermediate level (not existing elements)
     /            /           /           /            /
    /            /           /           /            /
   /            /           /           /            /
  --------------------------------------------------
          ---------------------------------------------------
         /     /      /     /     /     /     /     /      /
        /-----/------/-----/-----/-----/-----/-----/------/
       /     /      /     /     /     /     /     /      /
      /-----/------/-----/-----/-----/-----/-----/------/      fine level
     /     /      /     /     /     /     /     /      /
    /-----/------/-----/-----/-----/-----/-----/------/
   /     /      /     /     /     /     /     /      /
  --------------------------------------------------

```

### Example Rank Setup 1 ###

The first setup considers two particiapting ranks. Rank 0 holds the
coarse element and rank 1 all the fine elements.

```

          ---------------------------------------------------
         /     /      /     /  1  /                        /
        /-----/------/-----/-----/                        /
       /     /      /     /  1  /                        /
      /-----/------/-----/-----/           0            /
     /     /      /     /  1  /                        /
    /-----/------/-----/-----/                        /
   /     /      /     /  1  /                        /
  --------------------------------------------------

```

#### Layout on rank 0: ####

Please notice, that in this situation rank 0 is interpolating the
face information downwards to the fine level before sending the
face information to rank 1.
Overall, on rank 0 the sitation looks like this:

```

          ---------------------------------------------------
         /                        /                        /
        /                        /                        /
       /                        /                        /        coarse level
      /                        /                        /
     /                        /                        /
    /                        /                        /
   /                        /                        /
  --------------------------------------------------
     halo + from finer    /        fluid element
       element           /
                 from finer face

          ---------------------------------------------------
         /            /           /           /            /
        /            /           /           /            /
       /            /           /           /            /
      /------------/-----------/-----------/------------/   intermediate level
     /            /           /           /            /
    /            /           /           /            /
   /            /           /           /            /
  --------------------------------------------------
     halo + from finer    /       ghost from coarser
        elements         /            elements
                   from finer face

          ---------------------------------------------------
         /     /      /     /     /     /     /     /      /
        /-----/------/-----/-----/-----/-----/-----/------/
       /     /      /     /     /     /     /     /      /
      /-----/------/-----/-----/-----/-----/-----/------/           fine level
     /     /      /     /     /     /     /     /      /
    /-----/------/-----/-----/-----/-----/-----/------/
   /     /      /     /     /     /     /     /      /
  --------------------------------------------------
      halo elements       /      ghost from coarser
                         /          elements
                     send face

```

#### Layout on rank 1: ####

Rank 1 will receive the face information on the finest element level and
computes it.
However, it will NOT do any interpolation for the remote rank 0.

@note This is in contrast to the level descriptor where the elements
           are communicating on the coarsest element level. The faces will
           always communicate on the fineste level.

```

          -----------------------------
         /                        /
        /                        /
       /                        /                                 coarse level
      /                        /
     /                        /
    /                        /
   /                        /                        /
  --------------------------------------------------
     ghost from finer     /        not existing elem.
       element           /
                 no action for this face

          ---------------------------------------------------
         /            /           /           /            /
        /            /           /
       /            /           /
      /------------/-----------/--                          intermediate level
     /            /           /
    /            /           /
   /            /           /           /            /
  --------------------------------------------------
      ghost from finer    /       not existing elem.
        elements         /
                no action for this face

          ---------------------------------------------------
         /     /      /     /     /     /     /     /      /
        /-----/------/-----/-----/-----/-----/-----/------/
       /     /      /     /     /     /     /     /      /
      /-----/------/-----/-----/-----/-----/-----/------/           fine level
     /     /      /     /     /     /     /     /      /
    /-----/------/-----/-----/-----/-----/-----/------/
   /     /      /     /     /     /     /     /      /
  --------------------------------------------------
      fluid elements      /      halo from coarser
                         /          elements
             receive and compute face

```

### Example Rank Setup 2 ###

In this situation rank 0 holds the coarsest element.
The finer elements are distributed across two other ranks (1 and 2).

```

          ---------------------------------------------------
         /     /      /     /  2  /                        /
        /-----/------/-----/-----/                        /
       /     /      /     /  1  /                        /
      /-----/------/-----/-----/           0            /
     /     /      /     /  1  /                        /
    /-----/------/-----/-----/                        /
   /     /      /     /  1  /                        /
  --------------------------------------------------

```

#### Layout on rank 0: ####

Please notice, that in this situation rank 0 is interpolating the
face information downwards to the fine level before sending the
face information to rank 1 and 2.
Overall, on rank 0 the sitation looks like this:

```

          ---------------------------------------------------
         /                        /                        /
        /                        /                        /
       /                        /                        /        coarse level
      /                        /                        /
     /                        /                        /
    /                        /                        /
   /                        /                        /
  --------------------------------------------------
      ghost from finer    /        fluid element
       element           /
                 from finer face

          ---------------------------------------------------
         /            /           /           /            /
        /            /           /           /            /
       /            /           /           /            /
      /------------/-----------/-----------/------------/   intermediate level
     /            /           /           /            /
    /            /           /           /            /
   /            /           /           /            /
  --------------------------------------------------
      ghost from finer    /       ghost from coarser
        elements         /            elements
                   from finer face

          ---------------------------------------------------
         /     /      /     /     /     /     /     /      /
        /-----/------/-----/-----/-----/-----/-----/------/
       /     /      /     /     /     /     /     /      /
      /-----/------/-----/-----/-----/-----/-----/------/           fine level
     /     /      /     /     /     /     /     /      /
    /-----/------/-----/-----/-----/-----/-----/------/
   /     /      /     /     /     /     /     /      /
  --------------------------------------------------
      halo elements       /      ghost from coarser
                         /          elements
                     send face

```

#### Layout on rank 1: ####

Rank 1 will receive the face information on the finest element level and
computes it.
However, it will NOT do any interpolation for the remote rank 0.

@note This is in contrast to the level descriptor where the elements are
communicating on the coarsest element level. The faces will always communicate
on the finest level.

```

                               --------
                                  /
                                 /
                                /                                 coarse level
                               /
                              /
                             /
   /                        /                        /
  --------------------------------------------------
     not existing elem.   /   not existing elem.
                         /
                 no action for this face

          --                    ------
         /                        /
             not existing elem.  /
       /                        /
      /------------/-----------/--                          intermediate level
     /            /           /
    /            /           /
   /            /           /           /            /
  --------------------------------------------------
      ghost from finer    /       not existing elem.
        elements         /
                no action for this face

          --                    ------                   ----
                                  /                        /
        /-----/------/-----/-----/-----/-----/-----/------/
       /     /      /     /     /     /     /     /      /
      /-----/------/-----/-----/-----/-----/-----/------/           fine level
     /     /      /     /     /     /     /     /      /
    /-----/------/-----/-----/-----/-----/-----/------/
   /     /      /     /     /     /     /     /      /
  --------------------------------------------------
      fluid elements      /      halo from coarser
                         /          elements
             receive and compute face

```

#### Layout on rank 2: ####

Rank 2 will receive face information on the finest level, but is not doing
any interpolation for the remote rank 0.
The situation looks like this:

```

          ---------------------------------------------------
         /                        /                        /
        /                        /                        /
       /                        /                        /        coarse level
      /                        /                        /
     /                        /                        /
    /                        /                        /
   /                        /                        /
  --------------------------------------------------
    not existing elem.    /     not existing elem.
                         /
                 no action face

          ---        ----       ------       ----         ---
         /            /           /           /            /
                                 /
       /            /           /           /            /
      /--        --/-        --/--       --/--        --/   intermediate level
                  /           /           /            /
                             /
   /            /           /           /            /
  -----------------------------      ------      ---
      not existing elem.  /       not existing elem.
                         /
                   no action face

          ---------------------------------------------------
         /     /      /     /     /     /     /     /      /
        /-----/------/-----/-----/-----/-----/-----/------/
       /                        /                        /
      /---                  ---/--                    --/           fine level
     /                        /                        /
    /--                    --/--                    --/
   /                        /                        /
  --------------------------------------------------
      fluid element       /      halo from coarser
                         /          elements
                receive and compute face

```

More complex level jumps
------------------------

Of course, in general much more complex situations may arise in a realistic
mesh. However, these situations are covered by the upper analysis, as the
algorithm creates the face description recursively. So, from one level to
another level every situation should reduce to one of the upper cases. By
this approach any mesh situation can be handled.

