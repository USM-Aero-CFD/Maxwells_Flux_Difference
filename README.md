# Maxwells_Flux_Difference
This is a repository which solves the 2D and 3D Maxwell's equations using the RD flux-difference approach (method B). Another related numerical method provided here is the non-compact second-order finite-volume (method A) to benchmark the numerical results against the newly designed RD flux-difference method.

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
        ./Maxwells_Flux_Difference <option_1> <option_2> ... <option_11> <option_12>
    ```
2.  All the options to run this executable are listed in the table below:
    |    Options    |    Descriptions    |    Values    |
    |    ----    |    ----    |    ----    |
    |    1    |    to determine in this test case is in 2D or 3D    |    2D, 3D    |
    |    2    |    to determine if the default mesh will be generated    |    true, false    |
    |    3    |    input the grid number (is important if one is to generate the default mesh) *    |    integer    |
    |    4    |    choose if the generated mesh is to be randomized    |    true, false    |
    |    5    |    determine if the Maxwell's equations are to be computed in TM or TE mode    |    A (for TM), B (for TE)    |
    |    6    |    select the numerical method for computation    |    A (compact finite-volume), B (RD flux-difference)    |
    |    7    |    define the end time for the simulation    |    any float value    |
    |    8    |    define the time delta for the simulation    |    any float value (must be much smaller than the timeLast)    |
    |    9    |    input Gmsh file    |    Gmsh_2D.msh (for 2D) or Gmsh_3D.msh (for 3D)    |
    |    10    |    input file name for the element file **    |    element file (for 2D) or tetrahedron file (for 3D)    |
    |    11    |    input file name for the node vertex ***    |    node file name (for both 2D and 3D)    |
    |    12    |    cross section file name across y=0 for 3D    |    for 3D only    |

    - *    If one provides his / her own meshes, just an approximate value will be needed. However, it is still a _MUST_.
    - **   If one generates the default mesh using the authors' algorithms, the default would be _Element <elementNumber>.txt_ (2D) or _Tetrahedron <tetrahedronNumber>.txt_ (3D).
    - ***   If one generates the default mesh, the default node file name would be _Node <nodeNumber>.txt_.
3.  For example, if one wishes to generate the default mesh and compute using the `RD flux-difference` method, the following command should be invoked:
    ```bash
        ./Maxwells_Flux_Difference "2D" "true" "40" "false" "A" "B" "2.0" "0.01" "Gmsh_2D.msh" "Element 80.txt" "Node 80.txt"
    ```
    Its equivalence in _3D_ would be
    ```bash
        ./Maxwells_Flux_Difference "3D" "true" "40" "false" "A" "B" "1.0" "0.001" "Gmsh_3D.msh" "Tetrahedron 80.txt" "Node 80.txt" "CrossSection 80.txt"
    ```
4.  In short, one has to provide 9 parameters for 2D calculations, while 10 parameters for 3D calculations.

## Generate the default Meshes
1.  By setting the _<option_2>_ above as `true`, the executable file will run the `MeshConstruction*.h` class. This will generate the default mesh as depicted in the paper in the `build/` directory.
2.  If one wishes to generate the other meshes and uses them for computations, the mesh files should be placed in this `build/` folder, and the _<option_2>_ should be set as `false`. At the same time _<option_9>_, _<option_10>_ and _<option_11>_ in 2D case simulation (or _<option_9>_, _<option_10>_, _<option_11>_ and _<option_12>_ for 3D) shall be replaced accordingly.

## Citation
1.  For more numerical details, kindly refer to the paper **Compact second-order schemes for the electomagnetic waves on unstructured grids**.

## Results
## 2D Test Cases
1. Test case of second-order finite-volume (interpolation with computed gradient) on structured triangular mesh, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.05    |   0.01 |   -1.3010   |   -2.7268    |    0.4554    |
    |   0.0333  |   0.01 |   -1.4771   |   -3.0724    |    1.0537    |
    |   0.025  |   0.01 |   -1.6021   |   -3.3223    |    1.8826    |
2. Test case of second-order finite-volume (interpolation with computed gradient) on randomized grid, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.05    |   0.01 |   -1.3010   |   -1.6881    |    0.4661    |
    |   0.0333  |   0.01 |   -1.4771   |   -1.3104    |    1.0695    |
    |   0.025  |   0.01 |   -1.6021   |   -0.7652    |    1.8912    |
3. Test case of RD flux-difference on structured triangular mesh, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.05    |   0.01 |   -1.3010   |   -2.4089    |    0.1823    |
    |   0.0333  |   0.01 |   -1.4771   |   -2.7683    |    0.4322    |
    |   0.025  |   0.01 |   -1.6021   |   -3.0330    |    0.8177    |
4. Test case of RD flux-difference on randomized grid, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.05    |   0.01 |   -1.3010   |   -1.9692    |    0.1752    |
    |   0.0333  |   0.01 |   -1.4771   |   -2.2802    |    0.4310    |
    |   0.025  |   0.01 |   -1.6021   |   -2.5125    |    0.8167    |


## 3D Test Cases
1. Test case of second-order finite-volume (interpolation with computed gradient) on structured triangular mesh, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.10    |   0.01 |   -1.0000   |   -1.9805  |    21.2600    |
    |   0.0666  |   0.01 |   -1.1761   |   -2.2854   |    70.5290    |
    |   0.05  |   0.01 |   -1.3010   |    -2.5109   |    166.9060    |
2. Test case of second-order finite-volume (interpolation with computed gradient) on randomized grid, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.10    |   0.01 |   -1.0000   |   -1.5644  |    21.8625    |
    |   0.0666  |   0.01 |   -1.1761   |   -1.7404   |    72.7352    |
    |   0.05  |   0.01 |   -1.3010   |   -1.8318    |    166.1800    |
3. Test case of RD flux-difference on structured triangular mesh, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.10    |   0.01 |   -1.0000   |   -1.8429  |    3.3771    |
    |   0.0666  |   0.01 |   -1.1761   |   -2.1816   |    11.0352    |
    |   0.05  |   0.01 |   -1.3010   |   -2.4278    |    26.3233    |
4. Test case of RD flux-difference on randomized grid, with `timeLast = 0.5` using 1/3 Simpson's rule time-integration gives
    |   Grid Size   |   Time Delta   |   log10 deltaX   |   log10 L2   |    Execution time (s)    |
    |   ---    |    ---    |   ---   |    ---    |    ---    |
    |   0.10    |   0.01 |   -1.0000   |   -1.5959  |    3.4782    |
    |   0.0666  |   0.01 |   -1.1761   |   -1.8589   |    11.1394    |
    |   0.05  |   0.01 |   -1.3010   |   -2.0458    |    26.1008    |