APES Build INfrastructure
=========================

This repository collects the tools for building APES.
It is included as `bin` subrepository in most projects
and provides `bin/waf` for them.
There is also a `wscript` that provides some general
settings and projects can recurse into it.

Building waf
------------

To modify the waf file (for example after changing the
included tools), clone this repository into a directory
that is not within another project with a wscript and
a configured build. (waf will parse the directory
hierarchy to find the topmost wscript)

Make your changes (our additional scripts are found
in the `tools` subdirectory, third party stuff (like the
waf sources and fypp) is located in the `external`
subdirectory.
Then run `build_waf.sh` to generate a new waf file.
When committing changes to the sources please also update
the generated waf, as this is the easiest way to distribute
the updated file to all repositories.
