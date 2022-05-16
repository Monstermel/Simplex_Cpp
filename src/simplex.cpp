#include "simplex.h"

#include <Python.h>

#include <fstream>
#include <iostream>
#include <stack>

#include "config.h"
#include "variable.h"


namespace optimization {
    LinearProblem::LinearProblem(char const* _name)
        : name(_name), solution_dimension(0), changed_sign(false) {}

    LinearProblem::LinearProblem(const LinearProblem& linear_problem)
        : name(linear_problem.name),
          solution_dimension(linear_problem.solution_dimension),
          original_solution_dimension(linear_problem.original_solution_dimension),
          objective_function(linear_problem.objective_function),
          aux_objective_function(linear_problem.aux_objective_function),
          plt_objctv_fnctn(linear_problem.plt_objctv_fnctn),
          constraints(linear_problem.constraints),
          plt_cnstrnts(linear_problem.plt_cnstrnts),
          binary_constraints(linear_problem.binary_constraints),
          integer_constraints(linear_problem.integer_constraints),
          no_negative_constraints(linear_problem.no_negative_constraints),
          changed_sign(linear_problem.changed_sign),
          binary_problem(linear_problem.binary_problem),
          integer_problem(linear_problem.integer_problem),
          artificial_constrait(linear_problem.artificial_constrait),
          solution(linear_problem.solution),
          feasible(linear_problem.feasible),
          solution_value(linear_problem.solution_value) {
        std::vector<std::pair<size_t, Variable*>> idx;
        for (size_t i = 0; i < linear_problem.variables.size(); i++) {
            variables.push_back(linear_problem.variables[i]->clone());
            if (linear_problem.variables[i]->type == SPLITTED) {
                idx.push_back(
                    {i, static_cast<SplittedVariable*>(linear_problem.variables[i])->aux});
            } else if (linear_problem.variables[i]->type == AUXILIARY) {
                for (size_t j = 0; j < idx.size(); j++) {
                    if (idx[j].second == linear_problem.variables[i]) {
                        static_cast<SplittedVariable*>(variables[idx[j].first])->aux =
                            static_cast<AuxiliaryVariable*>(variables.back());
                        idx.erase(idx.begin() + j);
                        break;
                    }
                }
            }
        }
    }

    LinearProblem& LinearProblem::operator=(LinearProblem const& linear_problem) {
        name                        = linear_problem.name;
        solution_dimension          = linear_problem.solution_dimension;
        original_solution_dimension = linear_problem.original_solution_dimension;
        objective_function          = linear_problem.objective_function;
        aux_objective_function      = linear_problem.aux_objective_function;
        plt_objctv_fnctn            = linear_problem.plt_objctv_fnctn;
        constraints                 = linear_problem.constraints;
        no_negative_constraints     = linear_problem.no_negative_constraints;
        integer_constraints         = linear_problem.integer_constraints;
        binary_constraints          = linear_problem.binary_constraints;
        plt_cnstrnts                = linear_problem.plt_cnstrnts;
        changed_sign                = linear_problem.changed_sign;
        artificial_constrait        = linear_problem.artificial_constrait;
        integer_problem             = linear_problem.integer_problem;
        binary_problem              = linear_problem.binary_problem;
        solution                    = linear_problem.solution;
        feasible                    = linear_problem.feasible;
        solution_value              = linear_problem.solution_value;

        std::vector<Variable*>::iterator it;
        for (it = variables.begin(); it != variables.end(); it++) {
            if ((*it)->creator == this) {
                delete *it;
            }
        }
        variables.clear();
        std::vector<std::pair<size_t, Variable*>> idx;
        for (size_t i = 0; i < linear_problem.variables.size(); i++) {
            variables.push_back(linear_problem.variables[i]->clone());
            if (linear_problem.variables[i]->type == SPLITTED) {
                idx.push_back(
                    {i, static_cast<SplittedVariable*>(linear_problem.variables[i])->aux});
            } else if (linear_problem.variables[i]->type == AUXILIARY) {
                for (size_t j = 0; j < idx.size(); j++) {
                    if (idx[j].second == linear_problem.variables[i]) {
                        static_cast<SplittedVariable*>(variables[idx[j].first])->aux =
                            static_cast<AuxiliaryVariable*>(variables.back());
                        idx.erase(idx.begin() + j);
                        break;
                    }
                }
            }
        }
        return *this;
    }

    LinearProblem::~LinearProblem() {
        std::vector<Variable*>::iterator it;
        for (it = variables.begin(); it != variables.end(); it++) {
            if ((*it)->creator == this) {
                delete *it;
            }
        }
    }

    void LinearProblem::add_variable(Variable* variable) { variables.push_back(variable); }

    void LinearProblem::load_problem(std::ifstream& file) {
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

        // Procesamiento de variables no negativas
        size_t initial_solution_dimension = solution_dimension;
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


        if (binary_problem) {
            size_t aux = binary_constraints.size();
            for (size_t i = 0; i < aux; i++) {
                Mtrx eye = binary_constraints[i].coefficients;
                add_constraint(Constraint(eye, LESS_EQUAL, 1.0));
            }
        }

        // Creamos la funcion objetivo para la fase uno
        aux_objective_function      = objective_function;
        aux_objective_function.type = MIN;
        aux_objective_function.coefficients.setZero();
        // Transformamos las restricciones en igualdades
        for (it = constraints.begin(); it != constraints.end(); it++) {
            // Agreamos la variable de holgura o exceso
            std::vector<Constraint>::iterator mit;
            if (it->type != EQUAL) {
                for (mit = constraints.begin(); mit != constraints.end(); mit++) {
                    if (mit == it) {
                        switch (it->type) {
                            case LESS_EQUAL:
                                mit->add_column(1);
                                break;
                            case GREATER_EQUAL:
                                mit->add_column(-1);
                                break;
                            default:
                                break;
                        }
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
                aux_objective_function.add_column(0);
                // Agregamos la restriccion de no negatividad
                solution_dimension++;
                Mtrx eye                    = Mtrx::Zero(1, solution_dimension);
                eye(solution_dimension - 1) = 1;
                add_constraint(Constraint(eye, NOT_NEGATIVE, 0));
                // Agregamos la nueva variable
                std::stringstream variable_name;
                if (it->type == LESS_EQUAL) {
                    variable_name << "slack_";
                    variable_name << solution_dimension - original_solution_dimension;
                    variables.push_back(new SlackVariable(this, variable_name.str().c_str()));
                } else if (it->type == GREATER_EQUAL) {
                    variable_name << "excess_";
                    variable_name << solution_dimension - original_solution_dimension;
                    variables.push_back(new ExcessVariable(this, variable_name.str().c_str()));
                }
            }
            // Agregamos la variable artificial
            if (it->type == EQUAL || it->type == GREATER_EQUAL) {
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
                aux_objective_function.add_column(1);
                // Agregamos la restriccion de no negatividad
                solution_dimension++;
                Mtrx eye                    = Mtrx::Zero(1, solution_dimension);
                eye(solution_dimension - 1) = 1;
                add_constraint(Constraint(eye, NOT_NEGATIVE, 0));
                // Agregamos la nueva variable
                std::stringstream variable_name;
                variable_name << "artificial_";
                variable_name << solution_dimension - original_solution_dimension;
                variables.push_back(new ArtificialVariable(this, variable_name.str().c_str()));
            }
            // Actualizamos el tipo de restriccion
            it->type = EQUAL;
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

    ssize_t LinearProblem::is_optimal(Mtrx const& tableau, ObjectiveFunction_T type, bool ignore,
                                      ColumnSet const& _base) {
        ssize_t cof_idx = 0;
        for (ssize_t i = 0; i < tableau.cols() - 1; i++) {
            if (ignore && variables[i]->type == ARTIFICIAL && !_base.contains(i)) {
                continue;
            }
            if (type == MIN) {
                cof_idx = (tableau(0, cof_idx) < tableau(0, i)) ? i : cof_idx;
            } else {
                cof_idx = (tableau(0, cof_idx) > tableau(0, i)) ? i : cof_idx;
            }
        }
        if (type == MIN) {
            if (tableau(0, cof_idx) > 0.0) {
                if (tableau(0, cof_idx) < EPSILON) {
                    return -1;
                }
                return cof_idx;
            } else {
                return -1;
            }

        } else {
            if (tableau(0, cof_idx) < 0.0) {
                if (tableau(0, cof_idx) > -EPSILON) {
                    return -1;
                }
                return cof_idx;
            } else {
                return -1;
            }
        }
    }

    ssize_t LinearProblem::min_ratio(Mtrx const& tableau, size_t index) {
        feasible        = false;
        ssize_t min_idx = 0;  // Indice invalido temporalmente
        for (ssize_t i = 1; i < tableau.rows(); i++) {
            if (tableau(i, index) > 0) {
                if (!feasible) {
                    min_idx  = i;
                    feasible = true;
                } else {
                    min_idx = (tableau(i, tableau.cols() - 1) / tableau(i, index) <
                               tableau(min_idx, tableau.cols() - 1) / tableau(min_idx, index))
                                  ? i
                                  : min_idx;
                }
            }
        }
        return min_idx;
    }

    void LinearProblem::print_tableu_itr(Mtrx tableau, size_t itr) const {
        if (VERBOSE) std::cout << std::endl << "Tableau " << itr << ":\n" << tableau << std::endl;
        return;
    }

    bool LinearProblem::two_phase() {
        // creamos un problema temporal
        feasible                = false;
        LinearProblem tmp_prblm = *this;
        // Step 1
        for (std::vector<Constraint>::iterator it = tmp_prblm.constraints.begin();
             it != tmp_prblm.constraints.end(); it++) {
            if (it->value < 0.0) {
                it->coefficients *= -1.0;
                it->value *= -1.0;
                switch (it->type) {
                    case GREATER_EQUAL:
                        it->type = LESS_EQUAL;
                        break;
                    case LESS_EQUAL:
                        it->type = GREATER_EQUAL;
                        break;
                    default:
                        break;
                }
            }
        }
        // Step 2
        tmp_prblm.process_to_standard_form();
        // Step 4
        Mtrx tableau(tmp_prblm.constraints.size() + 1, tmp_prblm.solution_dimension + 1);
        // Funcion objetivo
        size_t n_clmns = tableau.cols();
        for (size_t i = 0; i < n_clmns; i++) {
            if (!(i + 1 < n_clmns)) {
                tableau(0, i) = 0.0;
            } else {
                tableau(0, i) = tmp_prblm.aux_objective_function.coefficients(i) * -1.0;
            }
        }
        // Restricciones
        size_t n_rows = tableau.rows();
        for (size_t i = 1; i < n_rows; i++) {
            for (size_t j = 0; j < n_clmns; j++) {
                if (!(j + 1 < n_clmns)) {
                    tableau(i, j) = tmp_prblm.constraints[i - 1].value;
                } else {
                    tableau(i, j) = tmp_prblm.constraints[i - 1].coefficients(j);
                }
            }
        }
        // Base inicial
        ColumnSet base;
        for (size_t i = 0; i < tmp_prblm.variables.size(); i++) {
            if (tmp_prblm.variables[i]->type == SLACK ||
                tmp_prblm.variables[i]->type == ARTIFICIAL) {
                base.insert(i);
            }
        }
        // Sacamos la variables artificiales de la funcion objetivo
        for (size_t i = 0; i < base.size(); i++) {
            if (tableau(0, base[i])) {
                long double aux = tableau(0, base[i]);
                for (size_t j = 1; j < n_rows; j++) {
                    if (tableau(j, base[i])) {
                        for (size_t k = 0; k < n_clmns; k++) {
                            tableau(0, k) = tableau(0, k) - aux * tableau(j, k);
                        }
                        break;
                    }
                }
            }
        }
        /*
        for (size_t i = 0; i < n_rows; i++) {
            long double mean = 0.0;
            for (size_t j = 0; j < n_clmns; j++) {
                mean += tableau(i, j);
            }
            mean /= tableau.cols();
            for (size_t j = 0; j < n_clmns; j++) {
                tableau(i, j) = tableau(i, j) / mean;
            }
        }
        */
        // Fase uno
        size_t  intrtn = 0;
        ssize_t enter_variable;
        while (print_tableu_itr(tableau, intrtn++),
               enter_variable = tmp_prblm.is_optimal(tableau, MIN, false, base),
               enter_variable != -1) {
            ssize_t left_variable = tmp_prblm.min_ratio(tableau, enter_variable);
            if (tmp_prblm.feasible) {
                // Actualizamos la base
                base[left_variable - 1] = enter_variable;
                // Hacemos una copia del tableau
                Mtrx old_tableau = tableau;
                //  Actualizamos el renglon de la variable que sale
                for (size_t i = 0; i < n_clmns; i++) {
                    tableau(left_variable, i) /= old_tableau(left_variable, enter_variable);
                }
                // Actualizamos el resto de la tabla
                for (ssize_t i = 0; i < tableau.rows(); i++) {
                    if (i == left_variable) {
                        continue;
                    }
                    for (size_t j = 0; j < n_clmns; j++) {
                        tableau(i, j) = old_tableau(i, j) -
                                        old_tableau(i, enter_variable) * tableau(left_variable, j);
                    }
                }
            } else {
                if (VERBOSE) puts("Problema infactible");
                return false;
            }
        }
        if (tableau(0, n_clmns - 1) > EPSILON) {
            if (VERBOSE) puts("Problema infactible");
            return false;
        }
        bool method_end = false;
        for (size_t i = 0; i < base.size(); i++) {
            if (tmp_prblm.variables[base[i]]->type == ARTIFICIAL) {
                method_end = true;
                break;
            }
        }
        if (!method_end) {
            // Fase dos, resolviendo sin quitar
            for (size_t i = 0; i < n_clmns; i++) {
                if (!(i + 1 < n_clmns)) {
                    tableau(0, i) = 0.0;
                } else {
                    tableau(0, i) = tmp_prblm.objective_function.coefficients(i) * -1.0;
                }
            }
            for (size_t i = 0; i < base.size(); i++) {
                if (tableau(0, base[i])) {
                    long double aux = tableau(0, base[i]);
                    for (size_t j = 1; j < n_rows; j++) {
                        if (tableau(j, base[i])) {
                            for (size_t k = 0; k < n_clmns; k++) {
                                tableau(0, k) = tableau(0, k) - aux * tableau(j, k);
                            }
                            break;
                        }
                    }
                }
            }
            while (
                print_tableu_itr(tableau, intrtn++),
                enter_variable = tmp_prblm.is_optimal(tableau, objective_function.type, true, base),
                enter_variable != -1) {
                ssize_t left_variable = tmp_prblm.min_ratio(tableau, enter_variable);
                if (tmp_prblm.feasible) {
                    // Actualizamos la base
                    base[left_variable - 1] = enter_variable;
                    // Hacemos una copia del tableau
                    Mtrx old_tableau = tableau;
                    //  Actualizamos el renglon de la variable que sale
                    for (size_t i = 0; i < n_clmns; i++) {
                        tableau(left_variable, i) /= old_tableau(left_variable, enter_variable);
                    }
                    // Actualizamos el resto de la tabla
                    for (ssize_t i = 0; i < tableau.rows(); i++) {
                        if (i == left_variable) {
                            continue;
                        }
                        for (size_t j = 0; j < n_clmns; j++) {
                            tableau(i, j) = old_tableau(i, j) - old_tableau(i, enter_variable) *
                                                                    tableau(left_variable, j);
                        }
                    }
                } else {
                    if (VERBOSE) puts("Problema infactible");
                    return false;
                }
            }
        }

        Mtrx pre_solution = Mtrx::Zero(1, tmp_prblm.solution_dimension);
        for (size_t i = 0; i < tmp_prblm.constraints.size(); i++) {
            pre_solution(base[i]) = tableau(i + 1, n_clmns - 1);
        }
        // Procesamos todas las variables
        tmp_prblm.solution                        = Mtrx::Zero(1, tmp_prblm.solution_dimension);
        std::vector<Variable*>::const_iterator it = tmp_prblm.variables.begin();
        for (size_t idx = 0; idx < tmp_prblm.solution_dimension && it != tmp_prblm.variables.end();
             it++, idx++) {
            (*it)->process(pre_solution, tmp_prblm.solution, idx);
        }

        tmp_prblm.solution_value = 0.0;
        for (size_t i = 0; i < solution_dimension; i++) {
            tmp_prblm.solution_value +=
                tmp_prblm.objective_function.coefficients(i) * tmp_prblm.solution(i);
        }
        // Cargamos los datos al problema original
        solution_value = tmp_prblm.solution_value;
        feasible       = tmp_prblm.feasible;
        solution       = Mtrx::Zero(1, tmp_prblm.original_solution_dimension);
        for (size_t i = 0; i < tmp_prblm.original_solution_dimension; i++) {
            solution(i) = tmp_prblm.solution(i);
        }
        return true;
    }

    void LinearProblem::branch_and_bound() {
        two_phase();
        log();
        if (feasible) {
            std::cout << "# Problema factible" << std::endl;
        } else {
            std::cout << "# Problema no factible" << std::endl;
        }
        print_solution();
        if (feasible) {
            // Para verificar que el programa no se rompa
            assert(static_cast<long long>(original_solution_dimension) ==
                   integer_constraints[0].coefficients.cols());
            // Unimos todas las restricciones de variable entera
            std::vector<Constraint>::iterator it;
            Mtrx                              int_vars = Mtrx::Zero(1, original_solution_dimension);
            for (it = integer_constraints.begin(); it != integer_constraints.end(); it++) {
                int_vars += it->coefficients;
            }
            // Stack para recorido LIFO
            std::stack<LinearProblem> tree;
            // Elegimos el pivote
            for (size_t i = 0; i < original_solution_dimension; i++) {
                // si la variable i debe ser entera y tiene solucion decimal
                if (int_vars(i) && fmodl(solution(i), 1.0)) {
                    Mtrx tmp_mtrx = Mtrx::Zero(1, original_solution_dimension);
                    tmp_mtrx(i)   = 1.0;
                    // Creamos el subproblema 2
                    LinearProblem tmp = *this;
                    tmp.add_constraint(
                        Constraint(tmp_mtrx, GREATER_EQUAL, floor(solution(i)) + 1.0));
                    tree.push(tmp);
                    // Creamos el subproblema 1
                    tmp = *this;
                    tmp.add_constraint(Constraint(tmp_mtrx, LESS_EQUAL, floor(solution(i))));
                    tree.push(tmp);
                    break;
                }
            }
            // Metodo de bnb
            long double bound{0};
            bool        bound_found = false;
            while (!tree.empty()) {
                LinearProblem tmp = tree.top();
                tree.pop();
                tmp.two_phase();

                tmp.log();
                if (feasible) {
                    std::cout << "# Problema factible" << std::endl;
                } else {
                    std::cout << "# Problema no factible" << std::endl;
                }
                tmp.print_solution();
                if (tmp.feasible) {
                    if (bound_found) {
                        if (objective_function.type == MAX && bound > tmp.solution_value) {
                            continue;
                        } else if (objective_function.type == MIN && bound < tmp.solution_value) {
                            continue;
                        }
                    }
                    bool factible_solution = true;
                    bool ignore;
                    for (size_t i = 0; i < original_solution_dimension; i++) {
                        if (int_vars(i) && fmodl(tmp.solution(i), 1.0)) {
                            factible_solution = false;
                            ignore            = false;
                            Mtrx tmp_mtrx     = Mtrx::Zero(1, original_solution_dimension);
                            tmp_mtrx(i)       = 1.0;
                            // Creamos el subproblema 2
                            LinearProblem tmp_0 = tmp;
                            Constraint    aux;
                            aux = Constraint(tmp_mtrx, GREATER_EQUAL, floor(tmp.solution(i)) + 1.0);
                            for (it = tmp.constraints.begin(); it != tmp.constraints.end(); it++) {
                                if ((*it) == aux) {
                                    factible_solution = true;
                                    ignore            = true;
                                    break;
                                }
                            }
                            if (ignore) continue;
                            tmp_0.add_constraint(aux);
                            tree.push(tmp_0);
                            // Creamos el subproblema 1
                            tmp_0 = tmp;
                            aux   = Constraint(tmp_mtrx, LESS_EQUAL, floor(tmp.solution(i)));
                            for (it = tmp.constraints.begin(); it != tmp.constraints.end(); it++) {
                                if ((*it) == aux) {
                                    factible_solution = true;
                                    ignore            = true;
                                    break;
                                }
                            }
                            if (ignore) continue;
                            tmp_0.add_constraint(aux);
                            tree.push(tmp_0);
                            break;
                        }
                    }
                    if (factible_solution) {
                        assert(solution.cols() == tmp.solution.cols());
                        if (!bound_found) {
                            solution_value = tmp.solution_value;
                            solution       = tmp.solution;
                            bound          = solution_value;
                            bound_found    = true;
                        } else {
                            if (objective_function.type == MAX && bound < tmp.solution_value) {
                                // Cota inferior
                                solution_value = tmp.solution_value;
                                solution       = tmp.solution;
                                bound          = solution_value;
                            } else if (objective_function.type == MIN &&
                                       bound > tmp.solution_value) {
                                // Cota superior
                                solution_value = tmp.solution_value;
                                solution       = tmp.solution;
                                bound          = solution_value;
                            }
                        }
                        // Solo se esta prendiendo la primera vez que encuentra una cota
                    }
                }
            }
        }
        return;
    }

    void LinearProblem::solve() {
        if (integer_problem) {
            branch_and_bound();
        } else {
            two_phase();
        }
        if (feasible) {
            print_solution();
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
        for (size_t i = 0; i < original_solution_dimension; i++) {
            std::cout << variables[i]->name << " =\t" << solution(i) << std::endl;
        }
        puts(
            "######################################################################################"
            "##############################################################################");
        /*
        if (!integer_problem && original_solution_dimension == 2) {
            plot();
        }
        */

        return;
    }

}  // namespace optimization
