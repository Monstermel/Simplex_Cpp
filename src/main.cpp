#include <iostream>

#include "simplex.h"

using namespace optimization;

int main(int argc, char* argv[]) {
    if (argc == 2) {
        Simplex problem("Simplex");
        try {
            problem.load_problem(argv[1]);
            problem.process_to_standard_form();
            problem.dual_simplex();
            problem.print_solution();
        } catch (char const* msg) {
            std::cout << "Error: " << msg << std::endl;
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } else {
        std::cout << "Error: No se proporciono un problema" << std::endl;
        return EXIT_FAILURE;
    }
}
