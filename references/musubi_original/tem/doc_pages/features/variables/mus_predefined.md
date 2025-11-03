title: Predefined Fortran functions in MUSUBI

On this page we would like to introduce a list of predefined fortran
functions that can be used in Musubi and a tutorial on how to use them.

# Predefined space time functions

To distinguish between spatial and transient function one uses the variable
table. Use `combined` to provide a temporal and a spatial table. Each can be
either a constant, a lua function, a predefined function
or a space time function. Here is an example:

```lua
variable = {
  {
    name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal = {
        ...
      }, -- temporal
      spatial = {
        ...
      } -- spatial
    } -- st_fun
  }, -- vel_x
} -- variable
```

# Predefined spatial funtions

* Parabolic profile
* Random distribution

## Parabolic profile

You can use the parabolic profile if you provide the shape of the parabol
and the amplitude.

Here is an example on how to you use it.

```lua
variable = {
  {
    name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal = {
        ...
      }, -- temporal
      spatial = {
        predefined ='parabol',
        shape = {
          kind = 'canoND',
          object = {
            center  = {-8.0, 0.0, 0.0},
            halfvec = {{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}
          } -- object
        }, -- shape
        amplitude = 1.0
      } -- spatial
    } -- st_fun
  }, -- vel_x
  ...
} -- variable
```

@todo: picture

## Random distribution

You can use random numbers for the spatial part.
Example:

```lua
variable = {
  {
    name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      temporal = {
        ...
      }, -- temporal
      spatial = {
        predefined ='random',
        min = 0.0, -- default
        max = 1.0 -- default
      } -- spatial
    } -- st_fun
  }, -- vel_x
  ...
} -- variable
```
@todo: picture

# Predefined temporal functions

* linear, smooth
* datafile
* cos

## Linear and Smooth

Here is an example on how to use the linear and smooth functionality for
spacetime functions.

```lua
variable = {
  {
    name = 'vel_x',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'combined',
      spatial = 1.0, -- as an example
      temporal = {
        predefined = {
          kind = 'linear', -- or 'smooth'
          minfactor = 0.0, -- default = 0.0
          maxfactor = 1.0, -- default = 1.0
          from_time = 0, -- default = 0.0
          to_time = 1000, -- default = 1.0
        } -- predefined
      } -- temporal
    } -- st_fun
  } -- vel_x
} -- variable
```
The following graphs show the differences between `smooth` and `linear`.

![linear_smooth](|media|/transient.png)

## Using data files

```lua
temporal = {
  kind = 'predefined',
  predefined = {
    kind = 'datafile',
    filename = '...',
    intp = '...',
    ramping = { -- possible
      rampVal = ..., -- ramping value
      rampT = ... -- ramping time
    },
    fac = 1.0, -- default
    periodic = 'false' -- default
  } -- predefined
} -- temporal
```

## Using cosine function

```lua
temporal = {
  kind = 'predefined',
  predefined = {
    kind = 'cos',
    frequency = 1.0, -- default
    phase = 0.0, -- default
    offset = 0.0 -- default
  } -- predefined
} -- temporal
```


