/*
    Cosas por agregar:
        Implementar la capicidad de que la funcion objetivo acepte constantes
*/

#include "simplex.h"

#include <cstdlib>
#include <fstream>
#include <regex>
#include <sstream>

#include "sys/stat.h"
#include "variable.h"

namespace optimization {


    Simplex::Simplex(char const* _name)
        : name(_name), solution_dimension(0), changed_sign(false), inverse_recalculation_rate(10) {}

    Simplex::~Simplex() {
        std::vector<Variable*>::iterator it;
        for (it = variables.begin(); it != variables.end(); it++) {
            if ((*it)->creator == this) {
                delete *it;
            }
        }
    }

    void Simplex::add_variable(Variable* variable) { variables.push_back(variable); }

    bool Simplex::has_solutions() const { return !overconstrained; }

    bool Simplex::is_unlimited() const { return unlimited; }

    bool Simplex::must_be_fixed() const { return has_to_be_fixed; }

    void Simplex::load_problem(char const* problem_name) {
        std::ifstream file(problem_name);

        if (file.is_open()) {
            // Variable enum para determiar que bloque estamos leyendo
            ParsingContext current_parsing_block = NONE;
            // Variable que cuenta el numero de variables registradas
            size_t current_var = 0;

            while (!file.eof()) {
                // Extraemos una linea hasta encontrar \n
                std::string buffer_init;
                getline(file, buffer_init);
                // Extraemos strings hasta encontrear un espacio
                // Si el primer caractar es un espacio lo extrae pero no lo guarda
                std::stringstream buffer(buffer_init);
                std::string       token;
                buffer >> token;

                if (token.length()) {
                    if (token == "Informacion:")
                        current_parsing_block = DATA;
                    else if (token == "Variables:")
                        current_parsing_block = VARS;
                    else if (token == "Restriciones:")
                        current_parsing_block = CONSTRAINTS;
                    else if (token == "Objetivo:")
                        current_parsing_block = OBJECTIVE;
                    else {
                        if (current_parsing_block == NONE) {
                            throw("Indentificador de bloque invalido");
                        }
                        switch (current_parsing_block) {
                            case DATA: {
                                try {
                                    solution_dimension = std::stoul(token);
                                } catch (...) {
                                    throw("Numero de variables invalido");
                                }
                            } break;

                            case VARS: {
                                // No podemos llegar aqui si DATA
                                if (!solution_dimension) {
                                    throw("No se ha cargado el numero de variables");
                                } else if (solution_dimension == current_var) {
                                    throw("Definicion inconsistente de variables");
                                }
                                Mtrx aux = Mtrx::Zero(1, solution_dimension);

                                aux(current_var) = 1;

                                std::string lower_bound;
                                lower_bound = token;
                                std::string variable_name;
                                buffer >> variable_name;
                                std::string upper_bound;
                                buffer >> upper_bound;

                                if (upper_bound.empty() || variable_name.empty() ||
                                    lower_bound.empty()) {
                                    throw("Definicion de variables invalida");
                                }

                                if (lower_bound != "inf") {
                                    long double aux_1;
                                    try {
                                        aux_1 = std::stold(lower_bound);
                                    } catch (...) {
                                        throw("Definicion de variables invalida");
                                    }
                                    if (aux_1 == 0.0) {
                                        add_constraint(Constraint(aux, NOT_NEGATIVE, 0.0));
                                    } else {
                                        add_constraint(Constraint(aux, GREATER_EQUAL, aux_1));
                                    }
                                }
                                if (upper_bound != "inf") {
                                    long double aux_1;
                                    try {
                                        aux_1 = std::stold(upper_bound);
                                    } catch (...) {
                                        throw("Definicion de variables invalida");
                                    }
                                    add_constraint(Constraint(aux, LESS_EQUAL, aux_1));
                                }

                                add_variable(new Variable(this, variable_name.c_str()));
                                current_var++;
                            } break;

                            case CONSTRAINTS: {
                                // No podemos llegar aqui sin DATA y VARS
                                if (!solution_dimension) {
                                    throw("No se ha cargado el numero de variables");
                                }
                                if (!current_var) {
                                    throw("No se han cargado las variables");
                                }
                                if (current_var != solution_dimension) {
                                    throw("Definicion inconsistente de variables");
                                }

                                Mtrx coefficients(1, solution_dimension);
                                try {
                                    coefficients(0) = std::stold(token);
                                } catch (...) {
                                    throw("Definicion de restricciones invalida");
                                }
                                for (size_t i = 1; i < solution_dimension; i++) {
                                    std::string aux;
                                    buffer >> aux;
                                    try {
                                        coefficients(i) = std::stold(aux);
                                    } catch (...) {
                                        throw("Definicion de restricciones invalida");
                                    }
                                }

                                std::string ct;
                                buffer >> ct;
                                std::string aux;
                                buffer >> aux;
                                if (aux.empty() || ct.empty()) {
                                    throw("Definicion de restricciones invalida");
                                }
                                long double bound;
                                try {
                                    bound = std::stold(aux);
                                } catch (...) {
                                    throw("Definicion de restricciones invalida");
                                }

                                if (ct == ">")
                                    add_constraint(Constraint(coefficients, GREATER_EQUAL, bound));
                                else if (ct == "<")
                                    add_constraint(Constraint(coefficients, LESS_EQUAL, bound));
                                else if (ct == "=")
                                    add_constraint(Constraint(coefficients, EQUAL, bound));
                                else {
                                    throw("Definicion de restricciones invalida");
                                }
                            } break;

                            case OBJECTIVE: {
                                if (solution_dimension == 0) {
                                    throw("No se ha cargado el numero de variables");
                                }
                                Mtrx costs(1, solution_dimension);
                                for (size_t i = 0; i < solution_dimension; i++) {
                                    std::string aux;
                                    buffer >> aux;
                                    try {
                                        costs(i) = std::stold(aux);
                                    } catch (...) {
                                        throw("Definicion de funcion objetivo invalida");
                                    }
                                }

                                if (token == "max")
                                    set_objective_function(ObjectiveFunction(MAX, costs));
                                else if (token == "min")
                                    set_objective_function(ObjectiveFunction(MIN, costs));
                                else {
                                    throw("Definicion de funcion objetivo invalida");
                                }
                            } break;
                            case NONE:
                                break;
                        }
                    }
                }
            }
            file.close();
            // Verificacion de lectura
            if (!solution_dimension || !current_var) {
                throw("Archivo de problema invalido");
            }
            if (current_var != solution_dimension) {
                throw("Archivo de problema invalido");
            }
            if (constraints.empty()) {
                throw("Archivo de problema invalido");
            }
            if (!(objective_function.coefficients.cols())) {
                throw("Archivo de problema invalido");
            }
            if (VERBOSE) {
                log();
            }
        } else {
            throw("No se pudo abrir el archivo");
        }

        return;
    }

    void Simplex::add_constraint(Constraint const& constraint) {
        if (solution_dimension != (size_t) constraint.coefficients.cols()) {
            throw("Error en el numero de las restricciones");
        }

        if (constraint.type == NOT_NEGATIVE) {
            nn_constraints.push_back(constraint);
        } else {
            constraints.push_back(constraint);
        }
    }

    void Simplex::set_objective_function(ObjectiveFunction const& objctv_fnctn) {
        if (solution_dimension != (size_t) objctv_fnctn.coefficients.cols()) {
            throw("Error en el tamaño de la funcion objetivo");
        }

        objective_function = objctv_fnctn;
    }

    void Simplex::log() const {
        std::cout << name << std::endl;
        std::cout << std::endl;

        std::cout << "Funcion objetivo:" << std::endl;
        objective_function.log();
        std::cout << std::endl;

        std::vector<Constraint>::const_iterator it;
        std::cout << "Restricciones: " << constraints.size() << std::endl;
        for (it = constraints.begin(); it != constraints.end(); it++) {
            it->log();
        }
        std::cout << std::endl;

        std::cout << "Restriciones de no negatividad: " << nn_constraints.size() << std::endl;
        for (it = nn_constraints.begin(); it != nn_constraints.end(); it++) {
            it->log();
        }
    }

    void Simplex::process_to_standard_form() {
        if (VERBOSE) {
            std::cout << std::endl << "# Transformando a forma estandar... ";
        }
        std::vector<Constraint>::iterator it;

        size_t initial_solution_dimension = solution_dimension;

        // Procesamiento de variables no negativas
        for (size_t i = 0; i < initial_solution_dimension; i++) {
            // Revisamos cada variable x_i
            bool has_constraint = false;

            // Determinamos si la variable x_i tiene restriccion de no negatividad
            for (it = nn_constraints.begin(); it != nn_constraints.end() && !has_constraint; it++) {
                if (it->coefficients(i)) {
                    has_constraint = true;
                }
            }

            // Si x_i es una variable libre la partimos en dos
            // x_i = x_i' - x_i''; x_i', x_i'' >= 0
            if (!has_constraint) {
                // Agregamos x_i'
                Mtrx eye = Mtrx::Zero(1, solution_dimension);
                eye(i)   = 1;
                add_constraint(Constraint(eye, NOT_NEGATIVE, 0));

                // Aumentamos el numero de variables
                solution_dimension++;
                std::vector<Constraint>::iterator mit;
                for (mit = nn_constraints.begin(); mit != nn_constraints.end(); mit++) {
                    mit->add_column(0);
                }

                // Agregamos x_i''
                Mtrx n_eye                    = Mtrx::Zero(1, solution_dimension);
                n_eye(solution_dimension - 1) = 1;
                add_constraint(Constraint(n_eye, NOT_NEGATIVE, 0));

                // Actualizamos las restricciones normales con la nueva variable x_i''
                for (mit = constraints.begin(); mit != constraints.end(); mit++) {
                    mit->add_column(-(mit->coefficients(i)));
                }

                // Actualizamos la funcion objetivo con la nueva variable x_i''
                objective_function.add_column(-(objective_function.coefficients(i)));

                // Creamos la variable separada y su auxiliar: x_i' y x_i''
                std::string aux_name(variables.at(i)->name);
                // Auxiliary sera la nueva variable osea x_i''
                Variable* auxiliary =
                    new AuxiliaryVariable(this, (aux_name + "_minus").c_str(), variables.size());
                // Splitted sera la nueva variable osea x_i' y tomara el lugar de x_i
                Variable* splitted = new SplittedVariable(this, variables.at(i)->name.c_str(),
                                                          (AuxiliaryVariable*) auxiliary);

                // Actualizamos las variables
                variables.at(i) = splitted;
                variables.push_back(auxiliary);
            }
        }

        // Procesamiento de las restricciones
        size_t aux = constraints.size();
        for (size_t i = 0; i < aux; i++) {
            if (constraints[i].type == EQUAL) {
                // Las restricciones del tipo
                // a_1x_1 + a_2x_2 + ... + a_nx_n = C
                // se transforman en dos
                // a_1x_1 + a_2x_2 + ... + a_nx_n <= C
                // a_1x_1 + a_2x_2 + ... + a_nx_n >= C

                // Cambiamos la restriccion de = a <=
                constraints[i].type = LESS_EQUAL;
                // Agregamos la restriccion auxiliar y la invertimos
                Mtrx eye = -1.0 * constraints[i].coefficients;
                add_constraint(Constraint(eye, LESS_EQUAL, -1.0 * constraints[i].value));
            } else if (constraints[i].type == GREATER_EQUAL) {
                // Invertimos las restricciones del tipo >=
                constraints[i].coefficients *= -1.0;
                constraints[i].value *= -1.0;
            }
        }
        for (it = constraints.begin(); it != constraints.end(); it++) {
            // Omitimos las restricciones ya procesadas
            if (it->type == EQUAL) {
                continue;
            }

            // Actualizamos el tipo de restriccion
            it->type = EQUAL;

            // Agreamos la variable de holgura
            std::vector<Constraint>::iterator mit;
            for (mit = constraints.begin(); mit != constraints.end(); mit++) {
                if (mit == it) {
                    mit->add_column(1);
                } else {
                    mit->add_column(0);
                }
            }
            for (mit = nn_constraints.begin(); mit != nn_constraints.end(); mit++) {
                mit->add_column(0);
            }

            // Ahora en la funcion objetivo
            objective_function.add_column(0);

            // Agregamos la restriccion de no negatividad
            solution_dimension++;
            Mtrx eye                    = Mtrx::Zero(1, solution_dimension);
            eye(solution_dimension - 1) = 1;
            add_constraint(Constraint(eye, NOT_NEGATIVE, 0));

            // Agregamos la nueva variable
            std::stringstream variable_name;
            variable_name << "slack_";
            variable_name << solution_dimension - 1;
            variables.push_back(new SlackVariable(this, variable_name.str().c_str()));
        }

        // Procesamiento de la funcion objetivo
        if (objective_function.type == MIN) {
            // Como usamos el metodo dual simplex y el simplex clasico
            // ocupamos siempre una funcion de maximizacion
            objective_function.type = MAX;
            changed_sign            = true;
            // Cambiamos el signo
            objective_function.coefficients *= -1.0;
        }

        // Comprobamos que todo se haya echo correctamente
        if ((size_t) objective_function.coefficients.cols() != solution_dimension) {
            throw("Error al transformar a la forma estandar");
        }
        if (variables.size() != solution_dimension) {
            throw("Error al transformar a la forma estandar");
        }
        for (it = constraints.begin(); it != constraints.end(); it++) {
            if ((size_t) (it->coefficients.cols()) != solution_dimension) {
                throw("Error al transformar a la forma estandar");
            }
        }
        for (it = nn_constraints.begin(); it != nn_constraints.end(); it++) {
            if ((size_t) (it->coefficients.cols()) != solution_dimension) {
                throw("Error al transformar a la forma estandar");
            }
        }

        // Actualizamos el nombre
        name += " (Forma estandar)";
        if (VERBOSE) {
            std::cout << "hecho" << std::endl << std::endl;
            log();
        }
    }
    /*
        void Simplex::solve() {
            ColumnSet initial_base;

            // Create alias to *this
            Simplex& original_problem = *this;

            // Create problem to work on
            Simplex standard_form_problem = original_problem;

            has_to_be_fixed = false;

            log();

            // Preprocessing
            if (VERBOSE) std::cout << "Generating problem in standard form ...";

            standard_form_problem.process_to_standard_form();

            if (VERBOSE) {
                std::cout << " done." << std::endl;
                standard_form_problem.log();
            }

            // Generate and solve artificial problem
            {
                // Create copy of standard form problem to create artificial problem
                Simplex artificial_problem = standard_form_problem;

                if (VERBOSE) std::cout << "Generating artificial problem ...";

                artificial_problem.process_to_artificial_problem();

                if (VERBOSE) {
                    std::cout << " done." << std::endl;
                    artificial_problem.log();
                }

                // Use artificial problem suggested base to solve it
                if (VERBOSE) std::cout << "Solving artificial problem ..." << std::endl;

                artificial_problem.solve_with_base(artificial_problem.suggested_base);

                if (VERBOSE) std::cout << "Done." << std::endl;

                if (artificial_problem.solution_value != 0) {
                    if (VERBOSE) std::cout << "Problem has no solution." << std::endl;

                    overconstrained = true;
                    return;
                } else {
                    overconstrained = false;

                    if (VERBOSE) {
                        std::cout << "Suggested initial base for original problem:";
                        artificial_problem.current_base.log(" ");
                    }

                    // If initial base doesn't contain artificial variables
                    // I can just use it, otherwise it may contain an artificial
                    // variable.

                    // Check for existence of a column index related  to an artificial
                    // variable by reading costs vector
                    int artificial_variable = -1;

                    for (int i = 0; i < artificial_problem.solution_dimension; ++i)
                        if (artificial_problem.objective_function.coefficients(i) == 1 &&
                            artificial_problem.current_base.contains(i))
                            artificial_variable = i;

                    // If index is still -1 (no artificial variables)
                    if (artificial_variable == -1) {
                        if (VERBOSE)
                            std::cout << "Base is clear about artificial variables, proceed ..."
                                      << std::endl;

                        standard_form_problem.suggested_base = artificial_problem.current_base;

                    } else {
                        //
                        //      If an artificial variable exists ... I can change the i (artificial)
                        //      column with a j column in current_out_of base so that//
                        //              *   j is not an auxiliary variable
                        //              *   (B^-1)_q * A^j != 0
                        //

                        if (VERBOSE)
                            std::cout << "Artificial variable detected in base: " <<
       artificial_variable
                                      << std::endl;
                        int  q = artificial_problem.current_base.index_of(artificial_variable);
                        Mtrx bi_row_q(1, (int) artificial_problem.current_base.size());

                        for (unsigned int k = 0; k < artificial_problem.current_base.size(); ++k)
                            bi_row_q(k) = artificial_problem.base_inverse(q, k);

                        // Find j
                        int j = -1;
                        for (unsigned int i = 0;
                             i < standard_form_problem.current_out_of_base.size() && j == -1; ++i) {
                            // Pick the ones that doesn't refer to an artificial variable
                            if (artificial_problem.costs(i) == 0) {
                                Mtrx column_j((int) standard_form_problem.current_base.size(), 1);

                                for (unsigned int k = 0; k <
       standard_form_problem.current_base.size();
                                     ++k)
                                    column_j(k) = artificial_problem.coefficients_matrix(k, i);

                                if ((double) (bi_row_q * column_j) != 0) j = i;
                            }
                        }

                        if (j != -1) {
                            // Found a j, substitute artificial_value with j
                            standard_form_problem.suggested_base = artificial_problem.current_base;
                            standard_form_problem.suggested_base.substitute(artificial_variable, j);
                            if (VERBOSE)
                                standard_form_problem.suggested_base.log("Now initial base is");

                        } else {
                            //    I didn't find a j which respected the requirements.
                            //    It may happen that for each j we have (B^-1)_q * A^j = 0,
                            //    this means that the rows of A are linearly dependent and
                            //    we can eliminate one of them. Let d be
                            //            d = e_q * B^-1
                            //    We have to eliminate a row for which d is non-zero.


                            std::cout << "Constraints are linearly dependent!" << std::endl;

                            // Find a constraint to eliminate (change)
                            int change = -1;
                            for (unsigned int i = 0;
                                 i < standard_form_problem.constraints.size() && change == -1; ++i)
                                if (bi_row_q(i) != 0) change = i;

                            std::cout << "Constraint #" << change << " must be eliminated."
                                      << std::endl;
                            has_to_be_fixed = true;
                            return;
                        }
                    }
                }
            }

            if (VERBOSE) std::cout << "Solving problem ..." << std::endl;

            standard_form_problem.solve_with_base(standard_form_problem.suggested_base);

            if (VERBOSE) std::cout << "Done." << std::endl;


            //      The solution of the standard form problem must be transposed to
            //      the original problem.

            if (VERBOSE) std::cout << "Processing standard form solution ..." << std::endl;

            if (standard_form_problem.unlimited) {
                unlimited = true;
            } else {
                unlimited = false;
                solution.resize(solution_dimension, 1);

                std::vector<Variable*>::const_iterator it;

                size_t index = 0;
                for (it = standard_form_problem.variables.begin();
                     it != standard_form_problem.variables.end(); it++, index++) {
                    (*it)->process(standard_form_problem.solution, solution, index);
                }

                solution_value     = standard_form_problem.solution_value;
                constraints_vector = standard_form_problem.constraints_vector;
                changed_sign       = standard_form_problem.changed_sign;
            }

            if (VERBOSE) std::cout << "Done." << std::endl;
        }
    */

}  // namespace optimization
