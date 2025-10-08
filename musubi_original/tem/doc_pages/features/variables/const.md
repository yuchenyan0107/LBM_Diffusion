title: Defining variables with constant values

A variable can be defined by a constant. This can be

* a Lua variable with constant value
* a Lua function result
* a constant value

@note Variables can be defined for [special shapes](../canonicalShapes.html). 
 
Here are examplex on how to define them:

```lua
function lua_fun1(x)
  value = math.cos(x)
  return value
end 

lua_fun2 = 1.0

variable = {
  {
    name = 'way1',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = 1.0
  },
  { 
    name = 'way2',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = {
      const = 1.0,
      shape = 'global'
    }
  },
  { 
    name = 'way3',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = lua_fun1(0)
  },
  {
    name = 'way4',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = lua_fun2
  }
}

-- or for a 3D vector

table.insert(variable,
  { name = 'way_3d',
    ncomponents = 3,
    vartype = 'st_fun',
    st_fun = {const = {0,0,0}}
  }
)
```
