#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

constexpr auto VERBOSE = false;

namespace optimization {
    typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> Mtrx;
    enum Constraint_T { LESS_EQUAL, GREATER_EQUAL, EQUAL, NOT_NEGATIVE };
    enum ObjectiveFunction_T { MAX, MIN };
    enum ParsingContext { NONE, DATA, VARS, CONSTRAINTS, OBJECTIVE };

    class Variable;
    class Constraint {
        friend class Simplex;

       public:
        Constraint(Mtrx const& cffcnts, Constraint_T _type, long double _value);

        void log() const;
        void add_column(long double _value);
        int  size() const;

       private:
        Mtrx         coefficients;
        Constraint_T type;
        long double  value;
    };
    class ColumnSet {
        friend class Simplex;

       public:
        void    insert(size_t column);
        void    remove(size_t column);
        void    log() const;
        bool    contains(size_t column) const;
        size_t  size() const;
        size_t& operator[](size_t idx);

       private:
        std::vector<size_t> columns;
    };
    class ObjectiveFunction {
        friend class Simplex;

       public:
        ObjectiveFunction();
        ObjectiveFunction(ObjectiveFunction_T _type, Mtrx const& cffcnts);
        ObjectiveFunction& operator=(ObjectiveFunction const& objective_function);

        // Solution value
        Mtrx get_value(Mtrx const& x) const;

        // Manipulation
        void add_column(long double value);

        // Debug
        void log() const;

       private:
        ObjectiveFunction_T type;
        Mtrx                coefficients;
    };
    class Simplex {
       public:
        Simplex(char const* _name);
        ~Simplex();

        void load_problem(char const* problem_name);
        void add_variable(Variable* var);
        void add_constraint(Constraint const& constraint);
        void set_objective_function(ObjectiveFunction const& objective_function);

        void log() const;
        void print_solution() const;
        bool has_solution() const;
        bool is_feasible() const;
        bool is_bounded() const;


        // Funcion que transforma a la forma estandar
        void process_to_standard_form();

        // El metodo implementado y funciones auxiliares
        void dual_simplex();

       protected:
        size_t is_optimal(Mtrx const& tableau);
        size_t min_ratio(Mtrx const& tableau, size_t index);

        // Datos del problema lineal
        std::string             name;
        size_t                  solution_dimension;
        ObjectiveFunction       objective_function;
        std::vector<Constraint> constraints;
        std::vector<Constraint> nn_constraints;
        std::vector<Variable*>  variables;
        bool                    changed_sign         = false;
        bool                    artificial_constrait = false;

        // Resultados
        Mtrx        solution;
        long double solution_value;
        bool        optimal = false;
        bool        feasible;
        bool        bounded;
    };

}  // namespace optimization

#endif
