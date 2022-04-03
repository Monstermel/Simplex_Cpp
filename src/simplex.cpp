/*
    Cosas por agregar:
        Implementar la capicidad de que la funcion objetivo acepte constantes
*/

#include "simplex.h"

#include <Python.h>

#include <cstdlib>
#include <ctgmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stack>

#include "sys/stat.h"
#include "variable.h"

namespace optimization {
    LinearProblem::LinearProblem(char const* _name)
        : name(_name), solution_dimension(0), changed_sign(false) {}

    LinearProblem::~LinearProblem() {
        std::vector<Variable*>::iterator it;
        for (it = variables.begin(); it != variables.end(); it++) {
            if ((*it)->creator == this) {
                delete *it;
            }
        }
    }

    void LinearProblem::add_variable(Variable* variable) { variables.push_back(variable); }

    void LinearProblem::load_problem(char const* problem_name) {
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
                    if (token == "Dimension")
                        current_parsing_block = DATA;
                    else if (token == "Variables")
                        current_parsing_block = VARS;
                    else if (token == "Restricciones")
                        current_parsing_block = CONSTRAINTS;
                    else if (token == "Objetivo")
                        current_parsing_block = OBJECTIVE;
                    else {
                        if (current_parsing_block == NONE) {
                            throw("Indentificador de bloque invalido");
                        }
                        switch (current_parsing_block) {
                            case DATA: {
                                try {
                                    solution_dimension          = std::stoul(token);
                                    original_solution_dimension = solution_dimension;
                                } catch (...) {
                                    throw("Definicion inconsistente de dimension");
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

                                std::string var_type;
                                var_type = token;
                                std::string variable_name;
                                buffer >> variable_name;

                                std::string lower_bound;
                                std::string upper_bound;
                                if (var_type != "bool") {
                                    buffer >> lower_bound;
                                    buffer >> upper_bound;
                                }

                                if ((upper_bound.empty() && var_type != "bool") ||
                                    (lower_bound.empty() && var_type != "bool") ||
                                    variable_name.empty() || var_type.empty()) {
                                    throw("Definicion de variables invalida");
                                }

                                if (var_type == "bool") {
                                    integer_problem = true;
                                    add_constraint(Constraint(aux, BINARY, 1.0));
                                    add_constraint(Constraint(aux, NOT_NEGATIVE, 0.0));
                                }

                                if (var_type == "int") {
                                    integer_problem = true;
                                    if (lower_bound == "0" && upper_bound == "1") {
                                        add_constraint(Constraint(aux, BINARY, 1.0));
                                        add_constraint(Constraint(aux, NOT_NEGATIVE, 0.0));
                                        var_type = "bool";
                                    } else {
                                        add_constraint(Constraint(aux, INTEGER, 0.0));
                                    }
                                }

                                if (var_type == "float") {
                                    binary_problem = false;
                                }

                                if (var_type != "bool" && lower_bound != "inf") {
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
                                if (var_type != "bool" && upper_bound != "inf") {
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

                                std::string cnstrnt_type;
                                buffer >> cnstrnt_type;
                                std::string cnstrnt_value;
                                buffer >> cnstrnt_value;
                                if (cnstrnt_value.empty() || cnstrnt_type.empty()) {
                                    throw("Definicion de restricciones invalida");
                                }
                                long double bound;
                                try {
                                    bound = std::stold(cnstrnt_value);
                                } catch (...) {
                                    throw("Definicion de restricciones invalida");
                                }

                                if (cnstrnt_type == ">")
                                    add_constraint(Constraint(coefficients, GREATER_EQUAL, bound));
                                else if (cnstrnt_type == "<")
                                    add_constraint(Constraint(coefficients, LESS_EQUAL, bound));
                                else if (cnstrnt_type == "=")
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

    void LinearProblem::add_constraint(Constraint const& constraint) {
        if (solution_dimension != (size_t) constraint.coefficients.cols()) {
            throw("Error en el numero de las restricciones");
        }

        if (constraint.type == NOT_NEGATIVE) {
            no_negative_constraints.push_back(constraint);
        } else if (constraint.type == INTEGER) {
            integer_constraints.push_back(constraint);
        } else if (constraint.type == BINARY) {
            binary_constraints.push_back(constraint);
        } else {
            constraints.push_back(constraint);
        }
    }

    void LinearProblem::set_objective_function(ObjectiveFunction const& objctv_fnctn) {
        if (solution_dimension != (size_t) objctv_fnctn.coefficients.cols()) {
            throw("Error en el tamaño de la funcion objetivo");
        }

        objective_function = objctv_fnctn;
    }

    void LinearProblem::log() const {
        std::cout << name << std::endl;
        std::cout << std::endl;

        std::cout << "Funcion objetivo:" << std::endl;
        objective_function.log();
        std::cout << std::endl;

        std::vector<Constraint>::const_iterator it;
        if (constraints.size()) {
            std::cout << "Restricciones: " << constraints.size() << std::endl;
            for (it = constraints.begin(); it != constraints.end(); it++) {
                it->log();
            }
            std::cout << std::endl;
        }

        if (no_negative_constraints.size()) {
            std::cout << "Restriciones de no negatividad: " << no_negative_constraints.size()
                      << std::endl;
            for (it = no_negative_constraints.begin(); it != no_negative_constraints.end(); it++) {
                it->log();
            }
            std::cout << std::endl;
        }

        if (integer_constraints.size()) {
            std::cout << "Variables enteras: " << integer_constraints.size() << std::endl;
            for (it = integer_constraints.begin(); it != integer_constraints.end(); it++) {
                it->log();
            }
            std::cout << std::endl;
        }

        if (binary_constraints.size()) {
            std::cout << "Variables binarias: " << binary_constraints.size() << std::endl;
            for (it = binary_constraints.begin(); it != binary_constraints.end(); it++) {
                it->log();
            }
            std::cout << std::endl;
        }
    }

    void LinearProblem::process_to_standard_form() {
        if (VERBOSE) {
            std::cout << "# Transformando a forma estandar... ";
        }
        // Creamos una copia para graficar
        if (solution_dimension == 2) {
            plt_cnstrnts     = constraints;
            plt_objctv_fnctn = objective_function;
        }

        std::vector<Constraint>::iterator it;

        size_t initial_solution_dimension = solution_dimension;

        // Procesamiento de variables no negativas
        for (size_t i = 0; i < initial_solution_dimension; i++) {
            // Revisamos cada variable x_i
            bool has_constraint = false;
            // Determinamos si la variable x_i tiene restriccion de no negatividad
            for (it = no_negative_constraints.begin();
                 it != no_negative_constraints.end() && !has_constraint; it++) {
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
                for (mit = no_negative_constraints.begin(); mit != no_negative_constraints.end();
                     mit++) {
                    mit->add_column(0);
                }
                for (mit = integer_constraints.begin(); mit != integer_constraints.end(); mit++) {
                    mit->add_column(0);
                }
                for (mit = binary_constraints.begin(); mit != binary_constraints.end(); mit++) {
                    mit->add_column(0);
                }

                // Agregamos x_i''
                Mtrx n_eye                    = Mtrx::Zero(1, solution_dimension);
                n_eye(solution_dimension - 1) = 1;
                add_constraint(Constraint(n_eye, NOT_NEGATIVE, 0));
                // Si ademas x_i es entera x_i'' tambien debe ser entera
                has_constraint = false;
                for (it = integer_constraints.begin();
                     it != integer_constraints.end() && !has_constraint; it++) {
                    if (it->coefficients(i)) {
                        has_constraint = true;
                    }
                }
                if (has_constraint) {
                    add_constraint(Constraint(n_eye, INTEGER, 0));
                }

                // Actualizamos las restricciones normales con la nueva variable
                // x_i''
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
        // if (!binary_problem) {
        //  Si el problema no es de tipo 0-1 entonces agregamos una restriccion del tipo
        //  x_k <= 1 donde x_k es una variable de tipo binaria
        aux = binary_constraints.size();
        for (size_t i = 0; i < aux; i++) {
            Mtrx eye = binary_constraints[i].coefficients;
            add_constraint(Constraint(eye, LESS_EQUAL, 1.0));
        }
        //}
        for (it = constraints.begin(); it != constraints.end(); it++) {
            // Omitimos las restricciones ya procesadas en iteraciones anteriores
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
            for (mit = no_negative_constraints.begin(); mit != no_negative_constraints.end();
                 mit++) {
                mit->add_column(0);
            }
            for (mit = integer_constraints.begin(); mit != integer_constraints.end(); mit++) {
                mit->add_column(0);
            }
            for (mit = binary_constraints.begin(); mit != binary_constraints.end(); mit++) {
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
            variable_name << solution_dimension - original_solution_dimension;
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
        aux = objective_function.coefficients.cols();
        for (size_t i = 0; i < aux; i++) {
            if (0 < objective_function.coefficients(i)) {
                // Agregamos una restriccion artificial
                artificial_constrait = true;
                // Agregamos su variable de holgura
                solution_dimension++;
                objective_function.add_column(0);
                std::vector<Constraint>::iterator mit;
                for (mit = constraints.begin(); mit != constraints.end(); mit++) {
                    mit->add_column(0);
                }
                for (mit = no_negative_constraints.begin(); mit != no_negative_constraints.end();
                     mit++) {
                    mit->add_column(0);
                }
                for (mit = integer_constraints.begin(); mit != integer_constraints.end(); mit++) {
                    mit->add_column(0);
                }
                for (mit = binary_constraints.begin(); mit != binary_constraints.end(); mit++) {
                    mit->add_column(0);
                }
                std::stringstream variable_name;
                variable_name << "slack_";
                variable_name << solution_dimension - original_solution_dimension;
                variables.push_back(new SlackVariable(this, variable_name.str().c_str()));
                // Restriccion de no negatividad
                Mtrx eye                    = Mtrx::Zero(1, solution_dimension);
                eye(solution_dimension - 1) = 1;
                add_constraint(Constraint(eye, NOT_NEGATIVE, 0));
                // Restriccion x_1 + ... + x_n + s_k = M
                size_t n_cnstrnts = constraints.size() + 1;
                for (size_t j = 0; j < solution_dimension - n_cnstrnts; j++) {
                    eye(j) = 1;
                }
                // Calculamos el valor de M
                long double buffer;
                long double max_value = 0.0;
                for (size_t j = 0; j < solution_dimension; j++) {
                    buffer    = std::abs(objective_function.coefficients(j));
                    max_value = (max_value < buffer) ? buffer : max_value;
                }
                for (size_t j = 0; j < constraints.size(); j++) {
                    for (size_t k = 0; k <= solution_dimension; k++) {
                        if (k == solution_dimension) {
                            buffer    = std::abs(constraints[j].value);
                            max_value = (max_value < buffer) ? buffer : max_value;
                        } else {
                            buffer    = std::abs(constraints[j].coefficients(k));
                            max_value = (max_value < buffer) ? buffer : max_value;
                        }
                    }
                }
                add_constraint(Constraint(eye, EQUAL, max_value * 100.0));
                break;
            }
        }

        // Comprobamos que todo se haya echo correctamente
        if ((size_t) objective_function.coefficients.cols() != solution_dimension) {
            throw("Error al transformar a la forma estandar 1");
        }
        if (variables.size() != solution_dimension) {
            throw("Error al transformar a la forma estandar 2");
        }
        for (it = constraints.begin(); it != constraints.end(); it++) {
            if ((size_t) (it->coefficients.cols()) != solution_dimension) {
                throw("Error al transformar a la forma estandar 3");
            }
        }
        for (it = no_negative_constraints.begin(); it != no_negative_constraints.end(); it++) {
            if ((size_t) (it->coefficients.cols()) != solution_dimension) {
                throw("Error al transformar a la forma estandar 4");
            }
        }

        // Actualizamos el nombre
        name += " (Forma estandar)";
        if (VERBOSE) {
            std::cout << "hecho" << std::endl << std::endl;
            log();
        }
    }

    size_t LinearProblem::is_optimal(Mtrx const& tableau) {
        size_t fixed_column = tableau.cols() - 1;
        size_t aux          = tableau.rows();
        size_t min_idx      = 1;
        for (size_t i = 2; i < aux; i++) {
            min_idx = (tableau(i, fixed_column) < tableau(min_idx, fixed_column)) ? i : min_idx;
        }
        if (tableau(min_idx, fixed_column) < 0.0) {
            return min_idx;
        } else {
            return 0;
        }
    }

    size_t LinearProblem::min_ratio(Mtrx const& tableau, size_t index) {
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

    bool LinearProblem::dual_simplex() {
        // Creacion del tableau
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
        // Hacemos una copia para evaluar el problema dual
        Mtrx dual_aux = tableau;

        // Creamos nuestra base inicial
        ColumnSet base;
        size_t    n_cnstrnts = constraints.size();
        for (size_t i = (n_clmns - 1) - n_cnstrnts; i < (n_clmns - 1); i++) {
            base.insert(i);
        }
        // Revisar si se añadio una restriccion artificial
        if (artificial_constrait) {
            size_t max_idx   = n_rows - 1;
            size_t max_idx_z = 0;
            for (size_t i = 1; i < n_clmns - 1; i++) {
                max_idx_z = (tableau(0, i) > tableau(0, max_idx_z)) ? i : max_idx_z;
            }
            // Actualizamos la base
            base[max_idx - 1] = max_idx_z;
            // Hacemos una copia del tableau
            Mtrx old_tableau = tableau;
            // Actualizamos el renglon de la variable que sale
            for (size_t i = 0; i < n_clmns; i++) {
                tableau(max_idx, i) /= old_tableau(max_idx, max_idx_z);
            }
            // Actualizamos el resto de la tabla
            for (size_t i = 0; i < n_rows; i++) {
                if (i == max_idx) {
                    continue;
                }
                for (size_t j = 0; j < n_clmns; j++) {
                    tableau(i, j) =
                        old_tableau(i, j) - old_tableau(i, max_idx_z) * tableau(max_idx, j);
                }
            }
        }

        // Metodo dual simplex
        std::cout << "\nInicio del metodo:";
        size_t intrtn = 0;
        size_t left_variable;
        while (std::cout << std::endl
                         << "Tableau " << intrtn << ":\n"
                         << tableau << std::endl,
               intrtn++, left_variable = is_optimal(tableau), left_variable) {
            size_t enter_variable = min_ratio(tableau, left_variable);
            if (feasible) {
                // Actualizamos la base
                base[left_variable - 1] = enter_variable;
                // Hacemos una copia del tableau
                Mtrx old_tableau = tableau;
                // Actualizamos el renglon de la variable que sale
                for (size_t i = 0; i < n_clmns; i++) {
                    tableau(left_variable, i) /= old_tableau(left_variable, enter_variable);
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
                puts("Problema infactible");
                return false;
            }
        }

        // El problema es factible o no acotado, asi que ahora evaluamos el dual
        size_t n_nbscs   = solution_dimension - n_cnstrnts;
        size_t d_n_rows  = n_nbscs + 1;
        size_t d_n_clmns = (artificial_constrait) ? n_cnstrnts + n_nbscs : n_cnstrnts + n_nbscs + 1;
        Mtrx   dual_tableau = Mtrx::Zero(d_n_rows, d_n_clmns);
        for (size_t i = 1; i < n_rows; i++) {
            if (artificial_constrait && !(i + 1 < n_rows)) {
                break;
            }
            dual_tableau(0, i - 1) = -dual_aux(i, n_clmns - 1);
        }
        for (size_t i = 0; i < n_nbscs; i++) {
            dual_tableau(i + 1, dual_tableau.cols() - 1) = -dual_aux(0, i);
        }
        for (size_t i = 1; i < n_rows; i++) {
            if (artificial_constrait && !(i + 1 < n_rows)) {
                break;
            }
            for (size_t j = 0; j < n_nbscs; j++) {
                dual_tableau(j + 1, i - 1) = -dual_aux(i, j);
            }
        }
        size_t aux_offset = (artificial_constrait) ? n_cnstrnts - 1 : n_cnstrnts;
        for (size_t i = 1; i < d_n_rows; i++) {
            dual_tableau(i, i + aux_offset - 1) = 1;
        }

        // Revisamos si necesita restriccion artificial
        long double max_value = dual_tableau(0, 0);
        for (size_t i = 0; i < d_n_clmns - d_n_rows; i++) {
            max_value = (max_value < dual_tableau(0, i)) ? dual_tableau(0, i) : max_value;
        }

        if (0.0 < max_value) {
            // Obtenemos el mayor valor en la tabla
            max_value = dual_tableau(0, 0);
            for (size_t i = 0; i < d_n_rows; i++) {
                for (size_t j = 0; j < d_n_clmns; j++) {
                    max_value = (max_value < dual_tableau(i, j)) ? dual_tableau(i, j) : max_value;
                }
            }

            // Agregamos la nueva restriccion
            Mtrx new_mtrx = Mtrx::Zero(d_n_rows + 1, d_n_clmns + 1);
            new_mtrx.topLeftCorner(d_n_rows, d_n_clmns - 1) =
                dual_tableau.topLeftCorner(d_n_rows, d_n_clmns - 1);
            new_mtrx.topRightCorner(d_n_rows, 1)  = dual_tableau.topRightCorner(d_n_rows, 1);
            d_n_rows                              = new_mtrx.rows();
            d_n_clmns                             = new_mtrx.cols();
            new_mtrx(d_n_rows - 1, d_n_clmns - 1) = max_value * 100;
            new_mtrx(d_n_rows - 1, d_n_clmns - 2) = 1;
            for (size_t i = 0; i < d_n_clmns - d_n_rows; i++) {
                new_mtrx(d_n_rows - 1, i) = 1;
            }
            dual_tableau = new_mtrx;

            // Forzamos la salida de la restriccion artificial
            size_t max_idx   = d_n_rows - 1;
            size_t max_idx_z = 0;
            for (size_t i = 1; i < d_n_clmns - 1; i++) {
                max_idx_z = (dual_tableau(0, i) > dual_tableau(0, max_idx_z)) ? i : max_idx_z;
            }

            // Hacemos una copia del tableau
            Mtrx old_tableau = dual_tableau;
            // Actualizamos el renglon de la variable que sale
            for (size_t i = 0; i < d_n_clmns; i++) {
                dual_tableau(max_idx, i) /= old_tableau(max_idx, max_idx_z);
            }
            // Actualizamos el resto de la tabla
            for (size_t i = 0; i < d_n_rows; i++) {
                if (i == max_idx) {
                    continue;
                }
                for (size_t j = 0; j < d_n_clmns; j++) {
                    dual_tableau(i, j) =
                        old_tableau(i, j) - old_tableau(i, max_idx_z) * dual_tableau(max_idx, j);
                }
            }
        }

        // Repetimos el metodo dual simplex
        while (left_variable = is_optimal(dual_tableau), left_variable) {
            size_t enter_variable = min_ratio(dual_tableau, left_variable);
            if (feasible) {
                // Hacemos una copia del tableau
                Mtrx old_tableau = dual_tableau;
                // Actualizamos el renglon de la variable que sale
                for (size_t i = 0; i < d_n_clmns; i++) {
                    dual_tableau(left_variable, i) /= old_tableau(left_variable, enter_variable);
                }
                // Actualizamos el resto de la tabla
                for (size_t i = 0; i < d_n_rows; i++) {
                    if (i == left_variable) {
                        continue;
                    }
                    for (size_t j = 0; j < d_n_clmns; j++) {
                        dual_tableau(i, j) = old_tableau(i, j) - old_tableau(i, enter_variable) *
                                                                     dual_tableau(left_variable, j);
                    }
                }
            } else {
                puts("Problema no acotado");
                return false;
            }
        }

        // Guardamos los resultados obtenidos
        if (changed_sign) {
            solution_value = tableau(0, n_clmns - 1);
        } else {
            solution_value = -tableau(0, n_clmns - 1);
        }
        Mtrx pre_solution = Mtrx::Zero(1, solution_dimension);
        for (size_t i = 0; i < n_cnstrnts; i++) {
            pre_solution(base[i]) = tableau(i + 1, n_clmns - 1);
        }
        // Retiramos todas las variables agregadas
        solution = Mtrx::Zero(1, solution_dimension);
        // Procesamos todas las variables
        std::vector<Variable*>::const_iterator it = variables.begin();
        for (size_t idx = 0; idx < solution_dimension && it != variables.end(); it++, idx++) {
            (*it)->process(pre_solution, solution, idx);
        }

        return true;
    }

    void LinearProblem::branch_and_bound() {}

    void LinearProblem::solve() {
        process_to_standard_form();
        if (integer_problem) {
            if (dual_simplex()) {
                std::cout << "Cota inicial: " << solution_value << std::endl;
                print_solution();
            }
        } else {
            // Aqui me gustaria agregar tambien el metodo de las dos fases
            if (dual_simplex()) {
                std::cout << "Solution_value: " << solution_value << std::endl;
                print_solution();
            }
        }
        return;
    }

    void LinearProblem::plot() const {
        // Variables auxiliares
        std::stringstream buffer;
        std::string       bffr_aux;
        std::cout << "Graficando...";

        // Creamos la configuracion inicial del script
        std::string script(
            "import numpy as np\n"
            "import matplotlib.pyplot as plt\n"
            "plt.axhline(color = 'black', linewidth = 0.93)\n"
            "plt.axvline(color = 'black', linewidth = 0.93)\n"
            "plt.xlabel('x')\n"
            "plt.ylabel('y')\n");

        // Configuramos el rango en el que vamos a graficar
        buffer.str(std::string());
        buffer << "x = np.linspace(" << solution(0) << " - 1000, " << solution(0) << " + 1000)\n";
        buffer << "plt.xlim(" << solution(0) << " - 1000, " << solution(0) << " + 1000)\n";
        bffr_aux = buffer.str();
        script += bffr_aux;

        // Generamos la funcion de cada restricccion
        size_t aux = plt_cnstrnts.size();
        for (size_t i = 0; i < aux; i++) {
            buffer.str(std::string());
            if (plt_cnstrnts[i].coefficients(1)) {
                buffer << "y" << i;
                buffer << " = ";
                buffer << "((" << plt_cnstrnts[i].value << ")-(" << plt_cnstrnts[i].coefficients(0)
                       << ")*x)/(" << plt_cnstrnts[i].coefficients(1) << ")\n";
                buffer << "plt.plot(x, y" << i;
            } else {
                buffer << "plt.axvline(x = "
                       << plt_cnstrnts[i].value / plt_cnstrnts[i].coefficients(0)
                       << ", color = next(plt.gca()._get_lines.prop_cycler)['color']";
            }
            buffer << ", label = '" << plt_cnstrnts[i].coefficients(0) << "x" << std::showpos
                   << plt_cnstrnts[i].coefficients(1) << std::noshowpos << "y";
            if (plt_cnstrnts[i].type == EQUAL) {
                buffer << "=";
            } else if (plt_cnstrnts[i].type == LESS_EQUAL) {
                buffer << "≤";
            } else {
                buffer << "≥";
            }
            buffer << plt_cnstrnts[i].value << "')\n";

            bffr_aux = buffer.str();
            script += bffr_aux;
        }

        // Generamos la funcion de la funcion objetivo
        buffer.str(std::string());
        buffer << "yz"
               << " = "
               << "((" << solution_value << ")-(" << plt_objctv_fnctn.coefficients(0) << ")*x)/("
               << plt_objctv_fnctn.coefficients(1) << ")\n"
               << "plt.plot(x, yz, label = 'Z=" << solution_value << "', color = 'black')\n";
        bffr_aux = buffer.str();
        script += bffr_aux;

        // Graficamos el punto de la solucion
        buffer.str(std::string());
        buffer << "plt.plot(" << solution(0) << ", " << solution(1) << ", 'o', label = '("
               << solution(0) << ", " << solution(1) << ")')\n";
        bffr_aux = buffer.str();
        script += bffr_aux;

        // Ejecutamos el script
        Py_Initialize();
        PyRun_SimpleString(script.c_str());
        PyRun_SimpleString(
            "plt.grid()\n"
            "plt.legend()\n"
            "plt.show()\n");
        std::cout << "hecho\n";
        Py_Finalize();
    }

    void LinearProblem::print_solution() const {
        std::cout << "Solucion:" << std::endl;
        std::cout << "Z  =\t" << solution_value << std::endl;
        for (size_t i = 0; i < solution_dimension; i++) {
            // if (variables[i]->type == ORDINARY || variables[i]->type == SPLITTED) {
            std::cout << variables[i]->name << " =\t" << solution(i) << std::endl;
            // }
        }
        if (solution_dimension == 2) {
            plot();
        }
        return;
    }

}  // namespace optimization
