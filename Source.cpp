# include <iostream>
# include <cstring>
# include "Compilation2D.h"
# include "Compilation3D.h"

int main (int argc, char** argv) {
    if (argc >= 10 && ((strcmp(argv[1], "2D") == 0) ? true : false) ) {

		/***********************************************/
		/*/
		argv[0] = Maxwell_Flux_Difference;	(function name)
        argv[1] = 2D;			        (string)
		argv[2] = Construct_Grid;			(boolean)
		argv[3] = Grid_Size;				(int)
		argv[4] = Randomization;			(boolean)
		argv[5] = TM_TE_Mode;				(char)
		argv[6] = Method;					(char)
        argv[7] = Gmsh_2D.msh;				(string)
		argv[8] = Element_File;				(string)
		argv[9] = Node_File;				(string)
		/**/
		/***********************************************/

		Compilation2D gridSizeL2;
		// gridSizeL2.gridIteration(false, 40, false, 'A', 'D', "Element 40.txt", "Node 40.txt");
		gridSizeL2.gridIteration( (strcmp(argv[2], "true") == 0) ? true : false,
                                    stod(argv[3]),
                                    (strcmp(argv[4], "true") == 0) ? true : false,
                                    argv[5][0],
                                    argv[6][0],
                                    argv[7], argv[8], argv[9]);

        return 0;
    } else if (argc >= 11 && ((strcmp(argv[1], "3D") == 0) ? true : false) ) {

		/***********************************************/
		/*/
		argv[0] = Maxwell_Flux_Difference;	(function name)
        argv[1] = 3D;			            (string)
		argv[2] = Construct_Grid;			(boolean)
		argv[3] = Grid_Size;				(int)
		argv[4] = Randomization;			(boolean)
		argv[5] = TM_TE_Mode;				(char)
		argv[6] = Method;					(char)
        argv[7] = Gmsh_3D.msh;				(string)
		argv[8] = Tetrahedron_File;			(string)
		argv[9] = Node_File;				(string)
		argv[10] = Cross_Section_File;		(string)
		/**/
		/***********************************************/

		Compilation3D gridSizeL2;
		// gridSizeL2.gridIteration(false, 40, false, 'A', 'D', "Tetrahedron 40.txt", "Node 40.txt", "CrossSection 40.txt");
		gridSizeL2.gridIteration( (strcmp(argv[2], "true") == 0) ? true : false,
                                    stod(argv[3]),
                                    (strcmp(argv[4], "true") == 0) ? true : false,
                                    argv[5][0],
                                    argv[6][0],
                                    argv[7], argv[8], argv[9], argv[10]);

        return 0;
    } else {
		cout << "Please provide enough information!";

        return 1;
    };
}
