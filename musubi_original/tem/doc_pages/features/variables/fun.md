title: Using Lua functions for space time functions

In this section we introduce the use on how to define 
variables with Lua functions.

Lua functions are defined inside the code like this:

```lua
function velocityX(x,y,z,t)
  velX = math.sin(x) * math.cos(y) * math.cos(z) * math.exp(-t/tD)
  return velX
end
```

To use them for the simulation run you need to add it to the variable table.

```lua
variable = {
  {
    name = 'lua_fun_1',
    ncomponents = 1
    vartype = 'st_fun',
    st_fun = velocityX
  },
  -- or in this way:
  {
    name = 'lua_fun_2',
    ncomponents = 1
    vartype = 'st_fun',
    st_fun = {
      fun = velocityX
    }
  },
  -- or for a special shape:
  {
    name = 'lua_fun_3',
    ncomponents = 1
    vartype = 'st_fun',
    st_fun = {
      fun = velocityX
      shape = { 
        -- example: line
        kind = 'canoND', 
        object = {
          origin = {0.0,0.0,0.0}, 
          vec = {1.0,0.0,0.0}, 
        } -- object
      } -- shape
    } -- st_fun
  },
}
```

@note For more possible shapes have a look at [Canonical Shapes](../canonicalShapes.html).

A combination of a spatial Lua function and a temporal Lua function is shown below.

```lua
a = 10.5

function x(t)
  y = math.cos(a*t)
  return y
end

function lua_fun(x,y,z)
  res = x*y + x*z*a
  return res
end

variable = {
  {
    name = 'two_lua',
    ncomponents = 1,
    vartype = st_fun,
    st_fun = {
      kind = 'combined',
      spatial = lua_fun,
      temporal = x
    }
  }
}
```
