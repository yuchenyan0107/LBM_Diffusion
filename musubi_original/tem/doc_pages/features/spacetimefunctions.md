title: Space-Time functions

# Space-time functions

Space-time functions allow the definition of arbitrary space-time dependent
functions, and the multiplication of a spatial function with a temporal function
as an often required special case.

Space-time functions might be predefined as Fortran functions, simple constants
or Lua functions.
They might return scalars or one-dimensional arrays.
If only a single value is expected from a Lua function, the Lua function is
supposed to return a scalar value. Irregardless if the function is invoked by an
array valued routine or a plain scalar one. Otherwise, the Lua function has to
return a table with the correct number of entries identified by position.

Note, that this makes the interface in the Lua script the same:

* whenever a single return value is expected, a plain scalar should be returned.
* otherwise a table should be returned


## Formal definition of space-time functions

Space-time functions can be written following some rules. Those rules are
noted using the [extended Backus-Naur form](https://en.wikipedia.org/wiki/Extended_Backus%E2%80%93Naur_form),
however, we skipped diving too deep into the definition tree.

```
stfun = "stfun = ", stfundef ;
stfundef = "{", singlestfundef, "}" | "{ ", "{", singlestfundef "}", ["{", singlestfundef' "}"], "}" | directdef ;
singlestfundef = detaileddef, {shapedef}
directdef = numerical_literal | luafunction | simplevectordef ;
simplevectordef = "{", numerical_literal, [ numerical_literal ], "}" ;
shapedef = "shape = { ... }" ;
detaileddef = scalardef | vectordef | luafunctiondef | predefineddef | stfunlistdef ;
scalardef = "const = ", numerical_literal ;
vectordef = "const = ", simplevectordef ;
luafunctiondef = "fun = " , luafunction ;
predefineddef = "predefined = { ", string_literal, { arguments }, "}" ;
  | "predefined = ", string_literal, ",", { arguments } ;
stfunlistdef = "multiples = {", stfundef, [ stfundef ], "link = ", [ first | add ], "}" ;
```

@todo Move the linked pages one level up, as they are not related to variables.

* combination of temporal and spatial function. Each can be a
    * [predefined Fortran function](variables/mus_predefined.html)
    * [constant](variables/const.html)
    * [Lua function](variables/fun.html)
* [constant](variables/const.html)
* [Lua function](variables/fun.html)
* [operation](variables/operation.html)

## Shape

A space time function can be defined for certain shapes of the simulated
domain. Default is global shape. See [Canonical Shapes](../canonicalShapes.html) for more information.
