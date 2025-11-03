-------------------------------------------------------------------------------
-- This is a sample configuration file to illustrate the general options
-- available via the treelm library.
-- It is supposed to be included in the dependent project configuration files
-- and should always work and represent the current status of the code.
-------------------------------------------------------------------------------

-- LOCAL PARAMETER
cube_length = 2.0
stop_restart = 1.0

-- The io_buffer_size specifies the amount of memory reserved for blocks of
-- temporary memory allocated for copies of data in MB.
-- (You do not have to specify it at all).
io_buffer_size = 1 -- default is 8 MB

-- Set a file or pseudo file to write undesired logging output to
-- (on each process, therefore it should be some local file or device).
-- Defaults to /dev/null if not provided, which should be reasonable
-- in most setups.
null_device = '/dev/null'

-- Communication pattern to use in parallel executions
commpattern = 'isend_irecv' -- default is simple isend-irecv exchange
--          = 'isend_irecv_overlap' allow overlapping of isends and irecvs
--          = 'typed_isend_irecv' use a MPI datatype to describe the data
--          = 'gathered_type' minimized number of blocks
--          = 'isend_irecv_mpimem' isend_irecv exchange with buffers
--                                 allocated by MPI
--          = 'overlap_mpimem' isend_irecv_overlap with buffers allocated by MPI

-- filename to dump performance measurement
timing_file = 'timing.res' -- default is 'timing.res'

-- Detailed timing information.
-- Will be written if the timer table is defined and the file name is not
-- an empty string.
-- Note: in contrast to the timing_file above, the file provided here
--       will be overwritten if it exists already!
timer = {
  file = 'timeinfo', -- Default, timing output is deactivated if set to ''

  -- For each timer you may define the level of detail, you would like to get
  -- printed for it. The default is to write a summary for each timer that shows
  -- the min, max and sum of the timer across all processes.
  --
  -- Provide each timer in a key-value pair in a table, where the first
  -- entry is the timer name (case-insensitive) and the second is the level of
  -- detail.
  details = {
    {'simloop', 'details'}, -- will be written to timeinfo_overall.details
    {'init', 'summary'},    -- this is the default
    {'preciceWrite', 'ignore'} -- this timer will not be written
  }
  -- Note that it does not result in an error if you provide timer names, that
  -- do not exist, but the timeinfo file will list them as unknown.
}

-- Print memory upon finalize?
-- This is obtained from the Linux pseudo file /proc/self/status.
-- It is only printed by the first process in a prallel run, and only provides
-- data on this first process.
printRuntimeInfo = true --default is true

-- During the definition of the neighborhood in treelm we make use of an
-- alltoall. However, we do not really need to communicate with all processes,
-- but rather only with our neighbors. (< ~ 30)
-- For large process counts it might get expensive to perform an alltoall in
-- this case. We implemented a sparse_alltoall, which communicates only with
-- the actual processes, with this flag you can activate or deactivate its
-- use. Typically, the sparse alltoall should not hurt for small process
-- counts and might be beneficial for large process counts.
--
-- Defaults to false.
use_sparse_alltoall = false

-- Settings for the primary logging object provided by treelm.
-- (Root outputs to stdout, if not file given)
-- The main control here is the level, which tells the logger how verbose
-- the output should be. The higher the log level, the more verbose your
-- output will be. Default is 1.
-- Use 0 to deactivate log messages in the primary logger.
-- Further settings are real_form to define the format of real numbers
-- (default = 'EN12.3') and int_form to define the format of integer
-- numbers (default = 'I0')
-- root_only controls whether only root should print the log messages.
-- If other processes shoul also write log messages, a non-empty filename
-- has to be provided, the processes will write to files with this name
-- and the process rank appended.
logging = {
  level = 1,
  filename = '',
  root_only = true,
  real_form = 'EN12.3',
  int_form  = 'I0',
}

--------------------------------------------------------------------------------
-- MESH DEFINITION
-- Defining the mesh can take two forms:
-- 1) simply a string specifying the directory to read the mesh from, if it is
--    indeed a directory there has to be a trailing "/", otherwise the string is
--    simply a prefix to all mesh files.
--mesh = 'mesh/'       -- (default = 'mesh/')
--weights = '' -- load weights for sparta load balancing
--             -- If this is an empty string '', no balancing will be done.
--             -- defaults to ''
--             -- Note: the actual file will get an endian suffix
--             --       appended (.lsb or .msb), depending on the
--             --       system.

-- File to write timer based weights to after the simulation.
-- Defaults to the same as weights, '' deactivates
-- writing of weights.
--write_weights = weights

--
-- 2) a table for a predefined mesh.
--    Available are fully periodic meshes:
--    - 'cube': a full cube, defined by origin, length and refinementLevel.
--    - 'slice': a plane out of the full cube.
--    - 'line': a line out of the full cube, if you use element_count instead
--              of refinementLevel for this mesh, an arbitrary number of
--              elements can be used for the mesh. In this case the length
--              refers to the overall length for the given number of elements
--              not the cube extent.
--    - 'line_bounded': a line with number of elements provided by
--                      element_count and boundaries in west and east.
--                      If this mesh is used, you need to provide a
--                      boundary condition for 'west' and 'east'.
mesh = {
  predefined = 'cube',      -- use the predefined full cube
  origin = {                -- origin of the cube
    (-1.0)*cube_length/2.0,
    (-1.0)*cube_length/2.0,
    (-1.0)*cube_length/2.0
  },
  length = cube_length,     -- length of the cube
  refinementLevel = 4       -- refinement level to resolve the cube
}
--------------------------------------------------------------------------------

-- Settings for load balancing 
-- To activate initial load balancing with weights based on level
-- balance = true -- activates simple initial load balancing 
-- or 
-- balance = { active = true,
--             balanceDir = 'balance/', -- directory where to store tmp data
--             kind = 'simple' -- based on predefined scalefactor for each level
--                             -- based on scaling (acoustic or diffusive)
--           --kind  = 'levelwise', -- weights based on level
--           --kind = 'boundary'-- weights computed from time taken to set boundary}
--             dynamic = true -- activates dynamic load balancing
--             time = {useIterations = true, -- Timings given in iterations? default: true
--                     max = tmax,           -- Maximal iteration to reach
--                     interval = tmax+1 }    -- Interval for checking density
--}
balance = false

--------------------------------------------------------------------------------
-- DEBUG OPTIONS
--debug = { 
--         debugMode = true,        -- default= false
--         verbose = 3,             -- default= 0

--         -- Use a logging object to track output in individual files.
--         -- The output messages to this logging object are usually
--         -- generated by tem_debug calls.
--         logging = { level = 1, -- how detailed should the log be?
--                     filename = 'dbg', -- to which file to write the log
--                     root_only = false }, -- should only root write msgs?

--         debugFiles = true,       -- default= false
--         -- What to dump into debugFiles? --
--           dumpTreeIDs = true,      -- default= false
--           dumpPropBits = true,     -- default= false
--           dumpAuxLists = true,     -- default= false
--           dumpDependencies = true, -- default= false
--           dumpState = true,        -- default= false
--           dumpHaloState = true,    -- default= false
--         --  end debugFiles  --
--         debugDependencies = true, -- default= false
--         checkDependencies = true, -- default= false
--         checkNans = true,         -- default= false
--         checkSteps = true,        -- default= false
--         debugMesh = 'dbg/mesh_',  -- default= ''
--         debugSource = true,       -- default= false
--         debugRestart = true,      -- default= false
--         traceMemory = true,       -- default= false
--        }
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- RESTART SETTINGS
restart = {
  -- Header of restart, to start simulation from.
  -- example:
  -- read = './restart/sample_lastHeader.lua',

  -- If an initialization should be done, when the file given in the
  -- read setting above is missing, set the init_on_missing to true.
  -- With an init_on_missing = false, the computation will stop, if
  -- the file in read is not found. This option has no meaning when no
  -- read is provided.
  -- default: false
  -- init_on_missing = false

  -- Prefix for the files, that are to be written during the
  -- simulation. If this ends with a path separator, the restart files
  -- will be written into the specified directory and that directory
  -- has to exist
  write = 'restart_', -- 'restart/' for a directory

  -- The points in time at which restart files are to be written is
  -- controlled by the time_control definition.
  -- Time definitions are provided in terms of:
  -- * sim:   simulation time, this is the time measurement for the
  --          transient phenomenon that is simulated.
  -- * iter:  the number of iterations that were done
  -- * clock: the real time that has passed since the beginning of
  --          the simulation in seconds.
  -- All three measures might be set, whatever is encountered first,
  -- will trigger the setting:
  -- {sim = 1.23, iter = 123, clock=12.3}
  -- If instead of a table a single number is provided, this is
  -- interpreted as the setting for sim.
  -- All not-given times are set to never.
  time_control = {
    -- Starting point after which restart files should
    -- be written.
    -- Setting iter to 0 here, results in restart files
    -- being written from the initial condition onwards.
    min = {iter = 0},

    -- The maximal point in time, up to which, restarts
    -- should be written.
    -- Note, that if this is not defined at all, it will
    -- be set to never, resulting in doing restarts for
    -- all times, including a final restart after
    -- reaching the termination of the time loop.
    max = stop_restart,

    -- Frequency at which restart files are to be
    -- written between min and max.
    interval = stop_restart,

    -- If a restart is read, it may be that the restart
    -- writing did not happen at the desired interval,
    -- nevertheless you might want to write the next
    -- restart file according to the original interval
    -- rythm. To achieve this you can tell the restart
    -- to align the trigger to last interval before
    -- the restart time with the align_trigger option:
    align_trigger = { sim = false },
    -- For each time component you can define whether
    -- an alignment should be done or not.
    -- If you only set one flag like:
    -- align_trigger = false
    -- It will be applied to the simtime, the others
    -- will be set to false.
    -- Default is false.

    -- Controlling if a restart file should be written
    -- or not involves communication, if a clock setting
    -- is used.
    -- This might have a negative impact on the
    -- performance, if it is done every iteration and
    -- a single iteration is too short.
    -- By setting check_iter the interval at which
    -- these checks are done, can be increased, and
    -- thus performance impacts reduced.
    -- CAUTION: be aware that increasing this to values
    --          larger than 1, decreases the accuracy
    --          of the points in time at which restarts
    --          are written!
    -- default:
    -- check_iter = 1
  }
}
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- TRACKING
--
-- tracking = {
--   {
--     label = 'label',                -- (default = 'unnamed')
--
--     variable = {'density', 'momentum'},
--
--     -- transient_reduction = {'sum','average'}, -- (default = 'none')
--
--     -- reduction = {'sum','average'},  -- (default = 'none')
--
--     shape = {
--       kind = 'canoND', 
--       object = {
--         origin = { centerX, centerY, centerZ }
--       }
--     },
--     -- See description in the restart settings for time_control.
--     time_control = {
--       min = {iter = 0},
--       max = stop_restart,
--       interval = {iter = 1}
--     }
--     folder = './tracking/'          --(default = './')
--     --for ateles if not ndof not specified it tries to dump the
--     --complete ndofs. Specify, ndofs = 1 for dumping the average
--     output = {
--       format = 'asciispatial',
--       ndofs = 1
--     }
--     -- It's also possible to dump the exact point value in ateles
--     -- for that you can set use_get_point to true, eg
--     -- output = {
--     __   format = 'ascii',
--     --   use_get_point = true
--     -- }
--     -- output = {format = 'harvester'}
--     -- output = {format = 'vtk'}
--   }
-- }

--------------------------------------------------------------------------------
-- TRANSIENT_REDUCTION
-- This is a reduction in time, so every value at every point in the tracking area is replaced in each timestep
-- by the result of a function that has probably access to data from previous timesteps.
-- This "transient-reduction" is applied before any spatial "reduction" operation.

-- REDUCTION
-- reduction = {'sum','average', 'l2norm', 'max', 'min', 'none'},
-- If we track more than one variable, there are two options to specify the reduction:
-- 1) give one reduction for all variables
-- 2) give as many reductions as there are are variables, the first reduction belongs to the first variable,
--    the second reduction to the second variable and so on.

-- SHAPE

-- Point
-- shape = {
--   kind = 'canoND', 
--   object = {
--     origin = {0.0, 0.0, 0.0} 
--   }
-- }

-- Line
-- shape = { 
--   kind = 'canoND', 
--   object = {
--     origin = {0, 0, 0}, 
--     vec = {0, 10, 0},
--     segments = {100}, 
--     distribution = 'equal' 
--   } 
-- }

-- Plane 
-- shape = {
--   kind = 'canoND',
--   object = {
--     origin = {0, 0, 0},
--     vec = {
--       {10, 0, 0}, 
--       {0, 10, 0}
--     }, 
--     segments = {100, 100}, 
--     distribution = 'equal' 
--   } 
-- }

-- Box 
-- shape = {
--   kind = 'canoND', 
--   object = {
--     origin = {0, 0, 0}, 
--     vec = {
--       {10, 0, 0}, 
--       {0, 10, 0},
--       {0, 0, 10}
--     }, 
--     segments = {100, 100, 100}, 
--     distribution = 'equal' 
--   } 
-- }
-- shape = { 
--   kind = 'canoND', 
--   object = {
--     origin = {0, 0, 0}, 
--     length = 10, 
--     segments = 100, 
--     distribution = ‘equal' 
--   } 
-- },

-- FORMAT
-- format = 'ascii',
-- format = 'asciiSpatial',
-- format = 'harvester',
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- POINT / LINE / PLANE DEFINITION

-- Point
-- shape = {
--   kind = 'canoND',
--   object = {
--     origin = {0.0, 0.0, 0.0} 
--   } 
-- }

-- Line 
-- shape = {
--   kind = 'canoND',
--   object = {
--     origin = {0, 0, 0}, 
--     vec = {0, 10, 0},
--     segments = {100}, 
--     distribution = 'equal' 
--   } 
-- }

-- Plane 
-- shape = { 
--   kind = 'canoND', 
--   object = {
--     origin = {0, 0, 0}, 
--     vec = {
--       {10, 0, 0}, 
--       {0, 10, 0}
--     }, 
--     segments = {100, 100}, 
--     distribution = 'equal' 
--   } 
-- }

-- Box 
-- shape = { 
--   kind = 'canoND', 
--   object = {
--     origin = {0, 0, 0}, 
--     vec = {
--       {10, 0, 0}, {0, 10, 0}, {0, 0, 10}
--     }, 
--     segments = {100, 100, 100}, 
--     distribution = 'equal' 
--   } 
-- }
-- shape = { 
--   kind = 'canoND', 
--   object = {
--     origin = {0, 0, 0}, 
--     length = 10, 
--     segments = 100, 
--     distribution = ‘equal' 
--   } 
-- },
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- PREDEFINED SPATIAL
---- gauss pulse
--{
--  predefined = 'gausspulse',
--  center = { 5.0, 5.0, 5.0 },    -- Gauss pulse center
--  halfwidth = 2.0,               -- half width of gauss pulse from center
--  amplitude = 2.0,               -- height or magnitude of gauss pulse
--  background = 1.225             -- reference value. In case of density, it is reference density
--]
--
---- 2d co-rotating vortex pair
--{
--  predefined = 'crvpX',                       -- or 'crvpY'
--  center = { 5.0, 5.0, 5.0 },    -- spinning center
--  radius_rot = 10.0,             -- distance of vortex centers / 2
--  circulation = 2.0              -- circulation of vortices
--]
--
---- parabol 2d
--{
--  predefined = 'parabol',
--  shape = { 
--    kind = 'canoND', 
--    object = {
--      origin = {0, 0, 0}, 
--      vec = {0, 10, 0},
--    }
--  }
--  amplitude = 1.0 -- {1.0,0.0,0.0} for vectorial variables
--}
--
---- parabol 3d
--{ 
--  predefined = 'parabol',
--  shape = { 
--    kind = 'canoND', 
--    object = {
--      origin = {0, 0, 0}, 
--      vec = { 
--        {10, 0, 0}, 
--        {0, 10, 0}
--      }, 
--    }
--  }
--  amplitude = 1.0 -- {1.0,0.0,0.0} for vectorial variables
--}
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
-- PREDEFINED TRANSIENT
---- smooth sin function
--{ 
--  predefined = 'smooth',
--  min_factor = 0.0,
--  max_factor = 1.0,
--  from_time = 0,
--  to_time = 1.0
--}
--
---- linear function
--{ 
--  predefined = 'linear',
--  min_factor = 0.0,
--  max_factor = 1.0,
--  from_time = 0,
--  to_time = 1.0
--}

----EXAMPLE
--
--velocity = { 
--  predefined = "combined",
--  spatial = {100,0.0,0.0} -- spatial const function
--  temporal = { 
--    predefined = 'smooth',  -- smooth sin function
--    min_factor = 0.0,
--    max_factor = 1.0,
--    from_time = 0,
--    to_time = 1.0
--  }
--}
--
---- in initial conditions:
--density = { 
--  predefined = 'gausspulse',
--  center = { 5.0, 5.0, 5.0 },    -- Gauss pulse center
--  halfwidth = 2.0,               -- half width of gauss pulse from center
--  amplitude = 2.0,               -- height or magnitude of gauss pulse
--  backgound = 1.225 }   -- reference value. In case of density, it is reference density
--}

--------------------------------------------------------------------------------
-- PREDEFINED SPACETIME

---- constant for all (x,y,z,t)
--spacetime = 5.0
--
---- combined, returns spatial(x,y,z)*transient(t)
--spacetime = {
--  spatial = { 
--    predefined = 'gausspulse',
--      center = { 5.0, 5.0, 5.0 },    -- Gauss pulse center
--      halfwidth = 2.0,               -- half width of gauss pulse from center
--      amplitude = 2.0,               -- height or magnitude of gauss pulse
--      backgound = 1.225
--    },
--  temporal = { 
--    predefined = 'smooth',  -- smooth sin function
--    min_factor = 0.0,
--    max_factor = 100.0,
--    from_time = 0,
--    to_time = 1.0
--  }
--}
--
---- predefined
--spacetime = { 
--  predefined = 'acoustic_pulse',
--  amplitude = 1.0,
--  halfwidth = 0.5,
--  speed_of_sound = 100.0,
--  background = 0.0,        -- (optional, default 0)
--  center = {0.0, 0.0, 0.0} -- (optional, default all 0)
--}
--
---- lua function
-- Space-Time functions can be expected to return scalars or 1d arrays.
-- Scalars are just single values that are returned, while arrays have
-- to be tables that are returned.
-- Note, that whenever there is only a single value expected, a scalar
-- should be returned, not a table with a single value!
--spacetime = luafunction
--------------------------------------------------------------------------------
