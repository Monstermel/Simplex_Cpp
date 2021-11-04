#include <Python.h>

#include <iostream>

#include "simplex.h"

using namespace optimization;

void plot(void) {
    Py_Initialize();
    PyRun_SimpleString("import numpy as np");
    PyRun_SimpleString("import matplotlib.pyplot as plt");
    PyRun_SimpleString("x = np.linspace(0, 15, 2000)");
    PyRun_SimpleString(
        "y1 = 0 * x + 2\n"
        "y2 = (25-x)/2.0\n"
        "y3 = (2*x-8)/4.0\n"
        "y4 = 2 * x -5\n"
        "yz = (73.75 - 4*x)/3");
    PyRun_SimpleString(
        "plt.plot(x, y1)\n"
        "plt.plot(x, y2)\n"
        "plt.plot(x, y3)\n"
        "plt.plot(x, y4)\n"
        "plt.plot(x, yz)\n"
        "y5 = np.minimum(y2, y4)\n"
        "y6 = np.maximum(y1, y3)\n"
        "plt.fill_between(x, y5, y6, where=y5>y6, color='grey', alpha=0.5)");
    PyRun_SimpleString("plt.show()");
    Py_Finalize();
    return;
}

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

        /*
        // Solve
        problem.solve();

        if (problem.must_be_fixed()) {
            std::cout << "El problema esta mal planteado" << std::endl;
            return EXIT_FAILURE;
        }

        if (problem.has_solutions()) {
            if (!problem.is_unlimited())
                problem.print_solution();
            else
                std::cout << "Problema con infinitas soluciones." << std::endl;
        } else {
            std::cout << "El problema es infactible" << std::endl;
        }
        */
        return EXIT_SUCCESS;
    } else {
        std::cout << "Error: No se proporciono un problema" << std::endl;
        return EXIT_FAILURE;
    }
}
