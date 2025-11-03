title: Operations inside the variable system

# Perform certain operations

It is possible to do operations with variables and store them in the
variable table in the prescribed way. It is defined in [[tem_variable_module]].

# What is available?

Here is a list of the currently available operations for the variable system.

* 'difference'

  Evaluate the function pointers of the dependent variables,
  and then calculate the difference between these two. ( scalar or vector )
  In lua file, first define new variable with varType operation `kind` as
  `difference` and provide two dependent variable via `input_varname`.
  If input_varname variable is not part of predefined solver variables then
  add also that variable via variable table for example spacetime function
  variable.
  For example: Define a variable called `dens_difference`, which depends on density
  and spacetime. One can get an error between simulation
  results and analytical solution.
  @note The number of input variables must be = 2.
         res(:) = a(:) - b(:)

```lua
variable = {
  {
    name = 'dens_reference',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = luaFun
  },
  {
    name = 'dens_difference',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'difference',
      input_varname = {'density', 'dens_reference'}
    }
  }
}
tracking = {
  variable = {'dens_difference'},
  folder = 'tracking/',
  shape = {
    kind = 'canoND',
    object = {
      origin = {3.0, 3.1, 3.0}
    }
  },
  format = 'ascii',
  time = {
    min = 0,
    max = tmax,
    interval = 1
  }
}
```

* 'rel_difference'

  @note Number of input variables = 2. res = ( a(:) - b(:) ) / b(:)

```lua
variable = {
  {
    name = 'dens_ref',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = luaFun
  },
  {
    name = 'dens_relDiff',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'rel_difference',
      input_varname = {'density', 'dens_ref'}
    }
  }
}
tracking = {
  variable = {'dens_relDiff'},
  folder = 'tracking/',
  shape = {
    kind = 'canoND',
    object = {
      origin = {3.0, 3.1, 3.0}
    }
  },
  format = 'ascii',
  time = {
    min = 0,
    max = tmax,
    interval = 1
  }
}
```

* 'addition'

  @note var3(:) = var1(:) + var2(:)

```lua
variable = {
  {
    name = 'var1',
    ncomponents = 3,
    vartype = "st_fun",
    st_fun = luaFun
  },
    name = 'var2',
    ncomponents = 3,
    vartype = "st_fun",
    st_fun = luaFun
  },
  {
    name = 'var3',
    ncomponents = 3,
    vartype = "operation",
    operation = {
      kind = 'addition',
      input_varname = {'var1', 'var2'}
    }
  }
}
```

* 'multiplication'

  Routine to multiply variables if all variables have same number of
  components.
  @note res(:) = a(:) * b(:)
         res(1) = a(1) * b(1)

```lua
variable = {
  {
    name = 'coeff',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = 0.25
  },
  {
    name = 'newVel',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'multiplication',
      input_varname = {'coeff', 'vel_mag'}
    }
  },
  ...
}
```

* 'division', 'div'

  Routine to divide variables if all variables have same number of
  components.
  @note res(:) = a(:) / b(:)

```lua
variable = {
  {
    name = 'coeff',
    ncomponents = 3,
    vartype = "st_fun",
    st_fun = 0.25
  },
  {
    name = 'newVel',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'division',
      input_varname = {'velocity', 'coeff'} -- numerator, denominator
    }
  },
  ...
}
```

* 'divide_vector_by_scalar'

  Routine to divide a vector by a scalar. So the second input variable must
  have one component.
  @note res(:) = a(:) / b

```lua
variable = {
  {
    name = 'coeff',
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = 0.25
  },
  {
    name = 'newVel',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'divide_vector_by_scalar',
      input_varname = {'velocity', 'coeff'} -- numerator, denominator
    }
  },
  ...
}
```

* 'gradient', 'grad'

  Calculates the gradient of the input variable. Only one input variable
  allowed. Number of components of the input variable can be 1,2 or 3.
  The resulting number of components must be ncomponents + 1.
  Only available in Ateles so far.

```lua
variable = {
  {
    name = 'dens_grad',
    ncomponents = 2, -- density: ncomponents=1
    vartype = "operation",
    operation = {
      kind = 'gradient',
      input_varname = 'density'
    }
  },
  ...
}
```

* 'magnitude'

  Evaluate magnitude of any vectorial variable. Number of components and
  number of input variables both have to be 1.
  In lua file, first define new variable with varType operation kind as
  "magnitude" and provide name of the variable from which magnitude
  to be derived in input_varname.
  If input_varname variable is not part of predefined solver variables then
  add also that variable via variable table.
  @note res = sqrt(sum( a(:) ))

```lua
variable = {
  {
    name = 'velMag',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'magnitude',
      input_varname = {'velocity'}
    }
  },
}
tracking = {
  variable = {'velMag'},
  folder = 'tracking/',
  shape = {
    kind = 'canoND',
    object = {
      origin = {3.0, 3.1, 3.0}
    }
  },
  format = 'ascii',
  time = {
    min = 0,
    max = tmax,
    interval = 1
  }
}
```

* 'extract'

  Extract component index of any vectorial variable.
  In lua file, first define new variable with varType operation kind as
  "extract" and provide name of the variable from which to extract
  component index via input_varname and
  index to extract via input_varIndex.
  If input_varname variable is not part of predefined solver variables then
  add also that variable via variable table.
  @note Both, ncomponents and number of input variables must be 1.

```lua
variable = {
  {
    name = 'vel_y',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'extract',
      input_varname = 'velocity',
      input_varindex = 2
    }
  }
}
```

* 'combine'

  Combine multiple variables into single variable with nComponent of
  output variable as sum of all input variables nComponents.
  In lua file, first define new variable with varType operation kind as
  "combine" and provide name of the variable from which to extract
  component index via input_varname (it must be single variable) and
  index to combine via input_varIndex.
  If input_varname variable is not part of predefined solver variables then
  add also that variable via variable table.

```lua
variable = {
  {
    name = 'dens_and_vel',
    ncomponents = 4,
    vartype = "operation",
    operation = {
      kind = 'combine',
      input_varname = {'density', 'velocity'}
    }
  },
}
```

* Boolean operationas and comparisons
  * 'greater_than', 'gt', '>'
  * 'greater_than_or_equal', 'ge', '>='
  * 'less_than', 'lt', '<'
  * 'less_than_or_equal', 'le', '<='
  * 'equal', 'eq', '='
  * 'not_equal', 'ne', '/='
  * 'and'
  * 'or'

  Using a boolean operation one will get either 1.0 (true) or 0.0 (false)
  as a resulting variable value. One must provide two input variables that
  must have the same number of components. Here is an example:

```lua
variable = {
  {
    name = 'var3',
    ncomponents = 3,
    vartype = "operation",
    operation = {
      kind = 'equal',
      input_varname = {'var1', 'var2'}
    }
  },
}

-- if var1 == var2 then:
-- var3 = {1.0, 1.0, 1.0}
```

* Temporal reduction [[tem_reduction_transient_module]]

  The temporal reduction operation allows the reduction of values over given
  iteration invervals. Thus it allows temporal averaging, as shown here:

```lua
  variable = {
    {
       name = 'press_timeavg',
       ncomponents = 1,
       vartype = "operation",
       operation = {
         kind='reduction_transient',
         input_varname={'pressure'},
         reduction_transient = {
           kind    = 'average',
           nrecord = 1000
         }
       }
    }
  }
```

  Other implemented reductions are sum, min and max.
