/*
    Fecha de entrega: 16/05/2022
    Cosas por hacer:
        Modificar el parser para que capture el objetivo del .problem
        Ademas de que pueda identificar entre problema tipo GRAPH o TRADITIONAL
        Implementar
            Min flow(con grafos)
*/

#include <fstream>
#include <iostream>

#include "graph.h"
#include "simplex.h"

using namespace optimization;

int main(int argc, char* argv[]) {
    (void) argc;
    std::ifstream file(argv[1]);
    if (file.is_open()) {
        std::string buffer_init;
        getline(file, buffer_init);
        std::stringstream buffer(buffer_init);
        std::string       token;
        buffer >> token;
        if (token == "GRAPH") {
            GraphProblem problem("Graph Problem");
            try {
                problem.load_problem(file);
                problem.solve();
            } catch (char const* msg) {
                std::cout << "Error: " << msg << std::endl;
                return EXIT_FAILURE;
            }
            return EXIT_SUCCESS;
        } else if (token == "SIMPLEX") {
            LinearProblem problem("Linear Problem");
            try {
                problem.load_problem(file);
                problem.solve();
            } catch (char const* msg) {
                std::cout << "Error: " << msg << std::endl;
                return EXIT_FAILURE;
            }
            return EXIT_SUCCESS;
        } else {
            std::cout << "Error: Algoritmo no seleccionado" << std::endl;
            return EXIT_FAILURE;
        }
    } else {
        throw("Error: No se pudo abrir el archivo");
        return EXIT_FAILURE;
    }
}
