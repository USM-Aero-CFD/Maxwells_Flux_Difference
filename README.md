# Maxwells_Flux_Difference
This is a repository which solves the 2D and 3D Maxwell's equations using the RD-flux-difference approach (method B). Other related numerical methods such as compact finite-volume (method A), RD-Galerkin (method C) and RD-LW (method D) are also given here.

## Build the project in CMake
1.  One has to install CMake (> Version 3.2) on his / her local computer and _git clone_ this project.
2.  Compiles the _C++_ project into a `build/` directory.
    ```bash
        cmake -S . -B build
    ```
3.  Builds the project into an executable file from the `build/` directory,
    ```bash
        cd build/ && cmake --build .
    ```

## Run the executable project
1.  Test the executable file from the `build/` directory with
    ```bash
        ./Maxwells_Flux_Difference <option_1> <option_2> ... <option_9> <option_10>
    ```
2.  All the options to run this executable are listed in the table below:
    |    Options    |    Descriptions    |    Values    |
    |    ----    |    ----    |    ----    |
    |    1    |    to determine in this test case is in 2D or 3D    |    2D, 3D    |
    |    2    |    to determine if the default mesh will be generated    |    true, false    |
    |    3    |    input the grid number (is important if one is to generate the default mesh) *    |    integer    |
    |    4    |    choose if the generated mesh is to be randomized    |    true, false    |
    |    5    |    determine if the Maxwell's equations are to be computed in TM or TE mode    |    A (for TM), B (for TE)    |
    |    6    |    select the numerical method for computation    |    A (compact finite-volume), B (RD flux-difference), C (RD-Galerkin), D (RD-LW)    |
    |    7    |    input Gmsh file    |    Gmsh_2D.msh (for 2D) or Gmsh_3D.msh (for 3D)    |
    |    8    |    input file name for the element file **    |    element file (for 2D) or tetrahedron file (for 3D)    |
    |    9    |    input file name for the node vertex ***    |    node file name (for both 2D and 3D)    |
    |    10    |    cross section file name across y=0 for 3D    |    for 3D only    |

    - *    If one provides his / her own meshes, just an approximate value will be needed. However, it is still a _MUST_.
    - **   If one generates the default mesh using the authors' algorithms, the default would be _Element <elementNumber>.txt_ (2D) or _Tetrahedron <tetrahedronNumber>.txt_ (3D).
    - ***   If one generates the default mesh, the default node file name would be _Node <nodeNumber>.txt_.
3.  For example, if one wishes to generate the default mesh and compute using the `RD flux-difference` method, the following command should be invoked:
    ```bash
        ./Maxwells_Flux_Difference "2D" "true" "40" "false" "A" "B" "Gmsh_2D.msh" "Element 40.txt" "Node 40.txt"
    ```
    Its equivalence in _3D_ would be
    ```bash
        ./Maxwells_Flux_Difference "3D" "true" "40" "false" "A" "B" "Gmsh_3D.msh" "Tetrahedron 40.txt" "Node 40.txt" "CrossSection 40.txt"
    ```
4.  In short, one has to provide 9 parameters for 2D calculations, while 10 parameters for 3D calculations.

## Generate the default Meshes
1.  By setting the _<option_2>_ above as `true`, the executable file will run the `MeshConstruction*.h` class. This will generate the default mesh as depicted in the paper in the `build/` directory.
2.  If one wishes to generate the other meshes and uses them for computations, the mesh files should be placed in this `build/` folder, and the _<option_2>_ should be set as `false`. At the same time _<option_7>_, _<option_8>_ and _<option_9>_ in 2D case simulation (or _<option_7>_, _<option_8>_, _<option_9>_ and _<option_10>_ for 3D) shall be replaced accordingly.

## Citation
1.  For more numerical details, kindly refer to the paper **Compact second-order schemes for the electomagnetic waves on unstructured grids**.
