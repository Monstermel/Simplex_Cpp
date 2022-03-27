#include <iostream>

#include "simplex.h"

using namespace optimization;

int main(int argc, char* argv[]) {
    if (argc == 2) {
        LinearProblem problem("Problem");
        try {
            problem.load_problem(argv[1]);
            problem.solve();
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
