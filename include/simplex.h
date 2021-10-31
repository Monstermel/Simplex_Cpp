#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <Eigen/Dense>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

constexpr auto TOL     = 0.00000000000000000001;
constexpr auto VERBOSE = true;

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
        void    substitute(size_t old_column, size_t new_column);
        int     index_of(size_t column) const;
        void    log(char const* prelude) const;
        bool    contains(size_t column) const;
        size_t& column(size_t idx);
        size_t  size() const;

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

        void solve();

        void print_solution() const;
        void log() const;

        bool is_unlimited() const;
        bool has_solutions() const;
        bool must_be_fixed() const;

        // Funcion que transforma a la forma estandar
        void process_to_standard_form();

       protected:
        // Column sets
        ColumnSet current_base;
        ColumnSet current_out_of_base;

        // Datos del problema lineal
        std::string             name;
        size_t                  solution_dimension;
        ObjectiveFunction       objective_function;
        std::vector<Constraint> constraints;
        std::vector<Constraint> nn_constraints;
        std::vector<Variable*>  variables;

        // Processed data
        Mtrx coefficients_matrix;
        Mtrx constraints_vector;
        Mtrx base_inverse;
        Mtrx column_p;

        size_t old_column;

        // Results
        Mtrx        base_solution;
        Mtrx        solution;
        Mtrx        reduced_cost;
        long double solution_value;

        bool optimal;
        bool unlimited;
        bool overconstrained;
        bool has_to_be_fixed;
        bool changed_sign;

        int inverse_recalculation_rate;
    };

}  // namespace optimization

#endif
