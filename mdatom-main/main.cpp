#include <MDSimulation.h>
#include <iostream>

#include <AbInitioInterface.h>

#include <cfenv>




/*!
 * main function.
 * The program expects one argument for the input file (parameters), one optional argument for
 * the file with the initial coordinates.
 */

int main(int argc, char *argv[]) {

    //Raises exception if there are invalid floating points e.g. 0/0 +-infinity sqrt(-1)
    feenableexcept(FE_INVALID | FE_OVERFLOW);

    switch (argc) {
        case 2:
            break;
        default:
            std::cerr << "Usage: mdatom input_file [coordinate_file] > output \n";
            return EXIT_FAILURE;
    }

    std::string parameterFile = argv[1];

    MDSimulation md(std::cout);
    try {
        md.performSimulation(parameterFile);
    }
    catch (std::exception &e) {
        std::cerr << e.what();
        return EXIT_FAILURE;
    }


    return EXIT_SUCCESS;
}
