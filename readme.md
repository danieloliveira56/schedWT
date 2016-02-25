Which cut separator to run is defined in file BB\BBNode.cpp

In order to use the OEC and Triangle Clique Cuts, the following folders should be created:

* OEC-Sep
* TriCliqueSep

Each folder content should be fetch from its corresponding git repository.

To compile, use:

* "make -f Makefile_Root" #compiles the executable swtLP-Root for solving only the Root Node
* "make -f Makefile_BCP" #compiles the executable swtLP-BCP for solving the problem to integrality or stop by time limit