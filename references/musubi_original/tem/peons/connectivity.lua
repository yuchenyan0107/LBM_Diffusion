-- How many neighbors should be searched for each element?
-- Ordering is according to the list of neighbors found in
-- tem_params_module. That is, the first 6 are the side
-- neighbors. The next 12 are the edge neighbors and the
-- final 8 are the corner neighbors.
neighbors = 6

-- File to write the resulting connectivity to in the form:
-- Element neighbor
-- with one line for each of these pairings.
outfile = 'connectivity.dat'
--------------------------------------------------------------------------------
-- MESH DEFINITION
cube_length = 2.0

-- Defining the mesh can take two forms:
-- 1) simply a string specifying the directory to read the mesh from, if it is
--    indeed a directory there has to be a trailing "/", otherwise the string is
--    simply a prefix to all mesh files.
--mesh = './mesh/'       -- (default = 'mesh/')
--
-- 2) a table for a predefined mesh. Right now only the full 'cube' without any
--    boundary conditions (periodic in all directions) is available as a
--    predefined mesh.
mesh = {
        predefined = 'cube',      -- use the predefined full cube
        origin = {                -- origin of the cube
                   (-1.0)*cube_length/2.0,
                   (-1.0)*cube_length/2.0,
                   (-1.0)*cube_length/2.0
                 },
        length = cube_length,     -- length of the cube
        refinementLevel = 1       -- refinement level to resolve the cube
       }
--------------------------------------------------------------------------------

-- Print memory upon finalize?
-- This is obtained from the Linux pseudo file /proc/self/status.
-- It is only printed by the first process in a prallel run, and only provides
-- data on this first process.
printRuntimeInfo = false --default is true

-- Settings for the primary logging object provided by treelm.
-- (Root outputs to stdout)
-- The main control here is the level, which tells the logger how verbose
-- the output should be. The higher the log level, the more verbose your
-- output will be. Default is 1.
-- Use 0 to deactivate log messages in the primary logger.
-- Further settings are real_form to define the format of real numbers
-- (default = 'EN12.3') and int_form to define the format of integer
-- numbers (default = 'I0')
logging = { level = 1 }
