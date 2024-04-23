# FormationControl

The following codes are used for computational simulations in the paper titled "Formation control between leader and migratory follower tissues allows coordinated growth of zebrafish midline tissues" by Kawanishi et al.


## Vertex Model for FP Simulation (for Mac)
This model is constructed using a programming code "main.cpp" and three header files in the folder named "Header". They are written in the C++ programming language. In our codes, we utilized the graphics library "GLSC3D" (https://github.com/GLSC3DProject/GLSC3D), a header file for calculating vectors (http://www3.u-toyama.ac.jp/akiyama/), and the Mersenne Twister random number generator (http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/emt19937ar.html) .　Due to the dependency of GLSC3D on the operating system, we have confirmed the following operations only on Mac.

First, install GLSC3D from GitHub. Folders named "bin", "include", and "lib" will be created under the Home directory. During the installation, you will be prompted to set the path in your shell, such as bash or z shell. Once the GLSC3D installation is complete, "glsc3d_3.h" will be stored in "include", and "libglsc3d_3.a" in "lib", preparing you for drawing graphics through GLSC3D. Next, download the Mersenne Twister from its website; you will obtain "mt19937ar.h" and "mt19937ar.c". Archive "mt19937ar.c" to create a file named "libmt19937ar.a". Save "mt19937ar.h" and "libmt19937ar.a" in the "include" and "lib" folders, respectively.

When you execute "makefile" in the folder via the command line, the executable file "run" will be generated. To execute "run", you need to specify three parameters: "alpha", "vmax", and "test_index", where the values for "alpha" and "vmax" are described in the Materials and Methods and Table S1, respectively, and "test_index" is an arbitrary integer. By default, we set alpha = 1, vmax = 2, and test_index = 1. Next, execute the executable file "run" by entering './run 1 2 1' into the command line.

## Vertex Model for Simulation with FP, HC, Notochord Cells (for Mac)
Similar to the above, this model is constructed using a programming code "main.cpp" and three header files in the folder named "Header". The execution procedure is the same as the one described above. Note that when the variable "Cell_Adhesion" in the "Parameter.h" file in the "Header" folder is uncommented, the posteriormost FP and HC cells are connected with the posteriormost notochord cells at specific vertices.
