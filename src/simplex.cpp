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


    Simplex::Simplex(char const* _name) : name(_name), solution_dimension(0), changed_sign(false) {}

    Simplex::~Simplex() {
        std::vector<Variable*>::iterator it;
        for (it = variables.begin(); it != variables.end(); it++) {
            if ((*it)->creator == this) {
                delete *it;
            }
        }
    }

    void Simplex::add_variable(Variable* variable) { variables.push_back(variable); }

    bool Simplex::has_solution() const { return optimal; }

    bool Simplex::is_feasible() const { return feasible; }

    bool Simplex::is_bounded() const { return bounded; }

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

    size_t Simplex::is_optimal(Mtrx const& tableau) {
        size_t fixed_column = tableau.cols() - 1;
        size_t aux          = tableau.rows();
        size_t min_idx      = 1;
        for (size_t i = 2; i < aux; i++) {
            min_idx = (tableau(i, fixed_column) < tableau(min_idx, fixed_column)) ? i : min_idx;
        }
        if (tableau(min_idx, fixed_column) < 0.0) {
            return min_idx;
        } else {
            optimal = true;
            return 0;
        }
    }

    size_t Simplex::min_ratio(Mtrx const& tableau, size_t index) {
        feasible       = false;
        size_t aux     = tableau.cols() - 1;
        size_t min_idx = 0;
        for (size_t i = 0; i < aux; i++) {
            if (tableau(index, i) < 0) {
                if (!feasible) {
                    min_idx  = i;
                    feasible = true;
                } else {
                    min_idx = (tableau(0, i) / tableau(index, i) <
                               tableau(0, min_idx) / tableau(index, min_idx))
                                  ? i
                                  : min_idx;
                }
            }
        }
        return min_idx;
    }

    void Simplex::dual_simplex() {
        Mtrx tableau(constraints.size() + 1, solution_dimension + 1);

        // Funcion objetivo y su valor
        size_t n_clmns = tableau.cols();
        for (size_t i = 0; i < n_clmns; i++) {
            if (!(i + 1 < n_clmns)) {
                tableau(0, i) = 0.0;
            } else {
                tableau(0, i) = objective_function.coefficients(i);
            }
        }
        // Restricciones y sus valores
        size_t n_rows = tableau.rows();
        for (size_t i = 1; i < n_rows; i++) {
            for (size_t j = 0; j < n_clmns; j++) {
                if (!(j + 1 < n_clmns)) {
                    tableau(i, j) = constraints[i - 1].value;
                } else {
                    tableau(i, j) = constraints[i - 1].coefficients(j);
                }
            }
        }
        // Creamos nuestra base inicial
        ColumnSet base;
        size_t    n_cnstrnts = constraints.size();
        for (size_t i = (n_clmns - 1) - n_cnstrnts; i < (n_clmns - 1); i++) {
            base.insert(i);
        }
        // Revisar si se añadio una restriccion artificial
        if (artificial_constrait) {
        }
        // Empezar metodo dual simplex
        size_t left_variable;
        while (left_variable = is_optimal(tableau), left_variable) {
            size_t enter_variable = min_ratio(tableau, left_variable);
            if (feasible) {
                Mtrx        old_tableau = tableau;
                long double pivot       = tableau(left_variable, enter_variable);
                // Actualizamos el renglon de la variable que sale
                for (size_t i = 0; i < n_clmns; i++) {
                    tableau(left_variable, i) /= pivot;
                }
                // Actualizamos el resto de la tabla
                for (size_t i = 0; i < n_rows; i++) {
                    if (i == left_variable) {
                        continue;
                    }
                    for (size_t j = 0; j < n_clmns; j++) {
                        tableau(i, j) = old_tableau(i, j) -
                                        old_tableau(i, enter_variable) * tableau(left_variable, j);
                    }
                }
            } else {
                throw("Problema infactible");
            }
        }
        std::cout << std::endl;
        std::cout << "Ultimo tableau:\n";
        std::cout << tableau;
        std::cout << std::endl;
        // Proceso los resultados obtenidos
    }

}  // namespace optimization
