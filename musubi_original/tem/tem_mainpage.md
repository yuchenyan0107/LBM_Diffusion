project: Treelm
summary: Tree based Elemental Mesh
project_website: https://apes-suite.github.io/pages/treelm.html
project_github: https://github.com/apes-suite/tem-source.git
src_dir: source
src_dir: utests
src_dir: build/ford
exclude: mpi_module.f90
exclude: tem_gamma_NAG.f90
exclude: tem_gamma_SLATEC.f90
exclude: tem_isNaN_vendor.f90
exclude: tem_isNaN_dummy.f90
external: aoturl = https://apes-suite.github.io/aotus
page_dir: doc_pages
media_dir: media
output_dir: docu
copy_subdir: media
graph: true
graph_maxdepth: 4
graph_maxnodes: 32
display: public
display: protected
display: private
sort: permission
source: true
externalize: true
print_creation_date: true
author: University of Siegen
md_extensions: markdown.extensions.toc

{!tem_index.md!}
