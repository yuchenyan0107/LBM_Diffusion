title: Tracking Module Description

@todo update

# Tracking

During the runtime of simulations, the user needs to keep track of the 
simultation state without being obliged to write out the complete domain 
to disk.
This is achieved by tracking objects, at which the state is evaluated 
at certain defined intervals and then written to disk.

## Tracking format

We support three formats for tracking.


- "ascii" 

Transient data: simple gnuplot-style 
A simple ascii file which creates a small header and then writes at each 
defined time step a column wise representation of all variables and points. 
This is not efficient and should only be used for simple point-tracking 
objects.
Each process writes its own file.

```lua
output = {
  format = 'ascii',
},
folder = './tracking/'
```

- "asciiSpatial" 

Spatial Data: This format stores solution data along with coordinates in 
single data file for every timestep. Therefore, there would be multiple 
files depending on the total number of time-steps. 
It is ideal for tracking *line* and *plane*.

```lua
output = {
  format = 'asciiSpatial',
},
folder = './tracking/'
```


- "Harvester Format"

Writes all participating elements to a treelm mesh and creates for each 
required time step a restart file. The output is done in `vtk` and `pvd` files
that can be visualised with `Paraview`.

It is also possible to configure the format for the timestamp to be
used in the filenames, see [[tem_timeformatter_load]].

This is a scalable tracking object which works regardless of number of 
processes.

~~~~~~~~~~~~~~~~~~~~~{.lua}
output = {
  format = 'vtk', 
},
folder = './output/'
~~~~~~~~~~~~~~~~~~~~~ 

@note We strongly encourage you to use the Harvester format. 
We will demonstrate each format in the examples of tracking point, 
line and plane in the following sections.

## Tracking variable types

- "variables" 
 
  All variables defined in the [[tem_varSys_module]] can be tracked.
  See [[tem_variable_module]] for how to define variables.
  @note Examples are given [here](variables/index.html).

- derived variables: all variables defined in the global variable system 
  can be tracked

- spacetime: as analytical result, can be defined in lua file

- [[tem_reduction_spatial_module]] reductions, e.g. difference, sum or
  average of defined variables

## Tracking shapes

We show the tracking with simple shapes here. 
Get an overview over all existing shapes in [[tem_shape_module]] and examples
are given [here](canonicalShapes.html)

- zero-dimensional: point 
- one-dimensional: line 
- two-dimensional: plane
- three-dimensional: cubes

So far point, line, plane are implemented in [[tem_shape_module]].

## Usage

The user specifies the desired parameters in the configuration file 
(in the scheme table) as

~~~~~~~~~~~~~~~~~~~~~{.lua}
-- if you want to have a analytical solution as reference, you can add 
--  a spacetime variable first
-- define single new (spacetime) variable
variable = { 
  {
    name='spacetime', 
    ncomponents=1, 
    kind = 'st_fun',
    st_fun = luaFun1 
  },
  {
    name='spacetime2', 
    ncomponents=1, 
    kind = 'st_fun',
    st_fun = luaFun2 
  },
  { 
    name='dens_diff', 
    ncomponents=1, 
    kind = 'operation',
    operation = {
      kind = 'difference'
      input_varname = {
        'density','spacetime'
      }
    }
  },
  { 
    name='wss_diff', 
    ncomponents=1, 
    kind = 'operation',
    operation = {
      kind = 'difference'
      input_varname = {
        'wss','spacetime2'
      }
    }
  }
}

-- Tracking
tracking = {
  { 
    label = 'probe_error',
    variable = {
      'velocity','dens_diff','wss_diff'
    },
    folder = './tracking/',
    output = {
      format = 'ascii'
    },
    time = {
      interval = 1, 
      min = 1000, 
      max = 4000 
    },
    reduction = {
      'l2norm',
      'l2norm',
      'l2norm'
    },  
    shape = { 
      kind = 'canoND',
      object = { 
        origin ={0.,0.,0.},
        vec = {1.,0.,0.}, 
        segments={100,100,100}
      }
    }
  }
} 
~~~~~~~~~~~~~~~~~~~~~

@note More examples for how to use the variable System can be found [here](variables/index.html).

### Reduction

The [[tem_reduction_spatial_module]] takes a set of input quantities and reduces them 
according to the chosen operation.
You can choose between the following operations:
- average: Average of all occurring values
- sum: Sum of all occurring values
- max: Maximum of all occurring values
- min: Minimum of all occurring values
- l2norm: L2 norm of all values, i.e. L_2(x) = sqrt( sum( x_i^2 ))
- linfnorm: L_inf norm of all values, i.e. L_inf(x) = max( | x_i | )

~~~~~~~~~~~~~~~~~~~~~{.lua}
reduction = {
              'norm' -- 'average' or 'sum'
            }
~~~~~~~~~~~~~~~~~~~~~

@note The table of reductions has to correspond to the number of entries in the
tracking variable table. For each variable, specify a reduction as

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = {
  variable = {'pressure', 'velMag', 'density'},
  reduction = { 'average', 'max', l2norm' }
}
~~~~~~~~~~~~~~~~~~~~~

### Using different shapes

Here, a point at the location {0.,0.,0.} is defined as a pressure probe. 
The interval is set to 1, meaning at each time step, the pressure has to 
be written.
The probing starts at simulation time 1000 and lasts until 4000.
The tracking write both pressure and velocityX in a single file.
Use ref_value to normalize the output.

```
# Name:       GaussPulse_TestSimulation     Date 2011-060-27  Time 11:43 CEST
# Probe type: Point at {0.0 0.0 0.0} 
# Quantities: pressure, velocityX
# tMin:       1000,     tMax:    4000,     Interval: 1
# ref_value: 1.0, 0.01
#   iTime     rTime      Pressure           velocityX
1000          1000.0     1.0241818          0.00158287
1001          1001.0     1.0251427          0.00152617
...                                      
4000          4000.0     1.7162511          0.01273673
``` 

The shape is defined as a canonical object. 
The different shapes are basically derivatives of a simple cube. 
A cube with no length in any direction is a point, a cube with length along 
only one axis is a line whereas a cube with length along two axes is simply 
treated as plane by the underlying code. 
One example each for Point, Line and Plane is shown below.

#### Point Tracking

<em>Example: tracking a point in ASCII format</em>

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = { label = 'Velocity', 
variable = {'velocity'},   -- options: density, velocity

--   Tracking over a point
shape={kind = 'canoND', origin = {0.5,0.5,0.5} },

time = {min = 0, interval = 1},     -- if no max given: active for all times
output = {format = 'ascii'},
folder = 'tracking/'}
}
~~~~~~~~~~~~~~~~~~~~~

Once, the above section is executed within the musubi.lua file, an output in 
the asciiTransient format is generated in the output folder specified 
(tracking in this particular example). 
Now, we open the result file inside tracking folder, the structure of which 
looks like below:

``` 
# Rank of the process:       0
#               time     velocity_01     velocity_02     velocity_03
1.0000  0.000000000000  0.000000000000  0.000000000000
2.0000  0.000000000000  0.000000000000  0.000000000000
3.0000  0.000000000000  0.000000000000  0.000000000000
4.0000  0.000000000000  0.000000000000  0.000000000000
...
...
1997.0000 -0.006065558329  0.000443400622 -0.000003481789
1998.0000 -0.006063901533  0.000443795703 -0.000003286949
1999.0000 -0.006061596302  0.000443327890 -0.000003534752
2000.0000 -0.006059665212  0.000443720604 -0.000003339771
```

Here, it can be noticed that 3 velocity components for all time steps at the 
point (0.5, 0.5, 0.5) are stored. The interval was chosen as 1, so they are 
tracked for each time step. If the interval was chosen as say 10, then the 
tracking would have been done after 10 units of simulated time.


#### Line Tracking

We need to set the origin and the direction vector, 
as e.g. {origin={0.,0.,0.}, vec={1.0,1.0,1.0}}.
The line is represented as a list of points in the tracking facility. 
The number of subdivisions and its distribution of the points on the line 
needs to be set, as e.g. `segments=500, distribution='equal'`. 
Make sure to specify enough segments to take into account all the required 
elements on the line. 
Elements are added only once, so it might be beneficial to specify much more 
segments than actually needed to consider elements on the highest levels as 
well.

<em>Example: Tracking a line in ASCII Spatial format</em>

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = { 
  label = 'line_probe',
  variable = {'state'}, -- state (=pdfs), pressure, density etc
  shape = { kind = 'canoND', 
            object={
              origin = {0.5,0.5,0.5}, 
              vec = { {0.0,0.2,0.0}},
              segments = {300},
              distribution='equal'} 
          },
  output = {format='asciispatial'}, 
  folder='tracking/',
  time = {interval = 5, min = 0, max = 800 } }
}
~~~~~~~~~~~~~~~~~~~~~

After, the above tracking is executed within musubi.lua file multiple res 
files are obtained inside the tracking folder each specific to particular 
time step. If you open any one of the 'res' files, there will be six data 
columns: three coordinate and three velocity components. 
An example is shown below: 

``` 
# Rank of the process:       0
#         coordX          coordY          coordZ     velocity_01     velocity_02     velocity_03
0.51562500      0.51562500      0.51562500 -0.000998448727 -0.000118182367  0.000007560959
0.51562500      0.54687500      0.51562500 -0.000922746216  0.000011023243 -0.000000130462
0.51562500      0.57812500      0.51562500 -0.001083956104  0.000038534363 -0.000003774649
...
0.51562500      0.64062500      0.51562500 -0.001343981219 -0.000131107478 -0.000004381864
0.51562500      0.67187500      0.51562500 -0.001261956611  0.000081027056  0.000000215364
0.51562500      0.70312500      0.51562500 -0.001733303068  0.000060154818  0.000005793207
``` 

Out of three coordinate components, only 'coordY' is varying and 'coordX', 
'coordZ' remain constant. This is according to our tracking input with 
`origin = {0.5, 0.5, 0.5}, vec = {0.0, 0.2, 0.0}`.
You can now make a line plot with any one of velocity components with 
respect to 'coordY' in [gnuplot](http://www.gnuplot.info/).

#### Plane Tracking

We need to set the origin and the direction vectors, as e.g. 
`object = {origin={0.,0.,0.}, vec={{1.0,1.0,1.0},{0.0,-1.0,1.0}}}`.
The plane is, just as the line, represented as a list of points in the 
tracking facility. 
The number of subdivisions in each direction and its distribution of the 
points on the plane needs to be set, 
as e.g. `segments1=500, segments2=300, distribution='equal'`. 
Make sure to specify enough segments to take into account all the required 
elements on the line. 
Elements are added only once, so it might be beneficial to specify much more 
segments than actually needed to consider elements on the highest levels 
as well.

<em> Example: tracking a plane in harvester format</em>

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = { 
  label = 'plane_probe',
  variable = {'state'}, -- state (=pdfs), pressure, density etc
  shape = { kind = 'canoND', 
            object={
              origin = {-75.,-15.,0.}, 
              vec = { {150.,0,0.}, {0.,30,0.}},
              segments = {300, 300},
              distribution='equal'} 
          },
  output = {format='vtk'}, 
  folder='tracking/',
  time = {interval = 5, min = 0, max = 800 } }
}
~~~~~~~~~~~~~~~~~~~~~ 

It will result in `vtk` files at the given time steps and one `pvd` file
that is an animation file over all `vtk` files.

#### Global mesh tracking

<em> Example: tracking the whole mesh</em>

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = { 
  ...
  shape = { kind = 'all' },
  -- or:
  shape = { kind = 'global' },
  -- or:
  shape = { kind = 'globalmesh' },
  ...
}
~~~~~~~~~~~~~~~~~~~~~

#### Boundary elements tracking

<em> Example: tracking the boundary elements</em>

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = { 
  ...
  shape = { kind = 'property', property = {'boundary'} },
  ...
}
-- or:

tracking = { 
  ...
  shape = { kind = 'boundary', boundary = {'bnd_name1','bnd_name2'} },
  ...
}
~~~~~~~~~~~~~~~~~~~~~

<em> Example: tracking solid elements</em>

~~~~~~~~~~~~~~~~~~~~~{.lua}
tracking = { 
  ...
  shape = { kind = 'property', property = {'solidified'} },
  ...
}
~~~~~~~~~~~~~~~~~~~~~

## Design and Implementation

The user-specified tracking shapes have to be translated to a set of
elements, on which the defined quantities are calculated and written to disk
in a suitable format.

All elements are converted to a set of points.
For each point, the correspoinding treeID is added to the element list of
the tracker as `track%elems( iLevel )%elem(1:nElems)`. 
Each element will only be added once.
All these elements are then treated in the main loop within the tracker
routine.

## peons/harvest_series.py

Each solver comes with a harvesting tool that allows the postprocessing
of restart data after the simulation has been run.
If this harvesting is to be done on a set of files, this may become a
little tedious.
To ease the post-processing there is a small Python script available in
`peons/harvest_series.py`, which takes care of the execution of the
post-processing script for multiple restart files.
It can be used as easily as running

~~~~~~~~~~~~~~~~~~~~~~~
peons/harvest_series.py restart_*.lua
~~~~~~~~~~~~~~~~~~~~~~~

The main argument this script expects is a list of files to work on.
There specification may make use of shell globbing patterns as in the
example above.

Aside from the files to process there are various other options that
can be specified, like the name of the harvesting executable, the
name of the Lua executable and the template to use to construct the
input file for the harvesting tool.
These options may always be provided as command-line arguments, but
can also be stored in a configuration file. By default harvest_series
will look for a file named series.config to read its configuration,
but you can tell it where to look for it via the `-c` option.

In the configuration file you provide the options to use in the form
of:

~~~~~~~~~~~~~~~~~~~~~~~
option: value
~~~~~~~~~~~~~~~~~~~~~~~

Thus, for example you could provide the name of the harvesting tool
as follows:

~~~~~~~~~~~~~~~~~~~~~~~
harvester: build/atl_harvesting
~~~~~~~~~~~~~~~~~~~~~~~

See the Python ConfigParser module documentation for more details on
the syntax of the configuration file.
You may provide a list of files to process in the config file. In this
case the files should be separated by commas, you may still use
globbing patterns in the file name specification.

Have a look at the available options with the `-h` option:

~~~~~~~~~~~~~~~~~~~~~~~
peons/harvest_series.py -h
~~~~~~~~~~~~~~~~~~~~~~~
