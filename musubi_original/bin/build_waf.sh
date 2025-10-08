set -x
UTESTS=$PWD/tools/utest_results.py
SEQF=$PWD/tools/command_sequence.py
FORD=$PWD/tools/make_fordoc.py
COCO=$PWD/tools/coco.py,$PWD/external/coco.f90
bindir=$PWD

TOOLS=clang_compilation_database,c_bgxlc,c_nec,doxygen,eclipse,fc_flang,fc_nag,fc_nec,fc_nfort,fc_bgxlf,fc_cray,fc_open64,fc_pgfortran,fc_solstudio,fc_xlf,print_commands,$UTESTS,$SEQF,$FORD
cd external/waf-2.1.5 && ./waf-light --make-waf --prelude='' --tools=$TOOLS,$COCO && cp waf $bindir && cd -
