TREe based ELemental Mesh
=========================

For an introduction to this library visit the
[Motivation page](page/motivation.html).
Treelm is released under the 2-clause BSD license please have a look into the
COPYRIGHT file for details.

The central definition of the treelmesh is found in the
[[treelmesh_module]].

Some details on the distributed handling of the tree are explained in the
[Description of the distributed octree structure](page/octree.html).

And the file format is described on the
[treelm file format page](page/fileformat.html).

Included Software
-----------------

Note, that this library includes some other software parts to your convenience.

- [Aotus](https://github.com/apes-suite/aotus) is used for configuration and
  header files. It is released under the terms of the MIT License as the
  [Lua library](http://www.lua.org), that it includes.
- [Waf](https://waf.io) is used as the build system and is released
  under the terms of the BSD license.
- [CoCo](http://users.erols.com/dnagle/coco.html) is used to preprocess Fortran
  source files and is released under the terms of the GPLv2.

License
=======

Treelm is licensed under the terms of the 2-clause BSD license reproduced below.
This means that Treelm is free software and can be used, reproduced, modified,
distributed and redistributed also for commercial purposes under the conditions
of the BSD license.
The only requirement is that some credit to the authors is given by putting this
copyright notice somewhere in your project.

According to good scientific practice, publications on results achieved in whole
or in part due to Treelm should cite at least one paper presenting the Treelm
software.

An appropriate reference is:

```
    @ARTICLE{Klimach:2012vi,
             author = {Klimach, Harald and Hasert, Manuel and Zudrop, Jens and Roller, Sabine},
             title = {{Distributed Octree Mesh Infrastructure for Flow Simulations}},
             journal = {European Congress on Computational Methods in Applied Sciences and Engineering},
             year = {2012},
             pages = {1--15}
    }
```
---

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY GERMAN RESEARCH SCHOOL FOR SIMULATION SCIENCES
GMBH, AACHEN “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GERMAN RESEARCH SCHOOL FOR
SIMULATION SCIENCES GMBH, AACHEN OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of German Research School for Simulation Sciences
GmbH, Aachen or University of Siegen.
