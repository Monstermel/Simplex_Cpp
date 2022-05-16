#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <Eigen/Dense>

namespace optimization {
    typedef Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic> Mtrx;
    enum Constraint_T { LESS_EQUAL, GREATER_EQUAL, EQUAL, NOT_NEGATIVE, INTEGER, BINARY };
    enum ObjectiveFunction_T { MAX, MIN };
    enum Variable_T { ORDINARY, SPLITTED, AUXILIARY, SLACK, EXCESS, ARTIFICIAL };
    enum ParsingContext { NONE, DATA, VARS, CONSTRAINTS, OBJECTIVE };

    class Variable;
    class Constraint {
        friend class LinearProblem;

       public:
        Constraint(){};
        Constraint(Mtrx const& cffcnts, Constraint_T _type, long double _value);

        bool operator==(const Constraint foo) {
            if (coefficients != foo.coefficients || type != foo.type || value != foo.value) {
                return false;
            }
            return true;
        }

        void log() const;
        void add_column(long double _value);
        int  size() const;

       private:
        Mtrx         coefficients;
        Constraint_T type;
        long double  value;
    };
    class ColumnSet {
        friend class LinearProblem;

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
        friend class LinearProblem;

       public:
        ObjectiveFunction();
        ObjectiveFunction(ObjectiveFunction_T _type, Mtrx const& cffcnts);
        ObjectiveFunction(const ObjectiveFunction& objective_function);
        ObjectiveFunction& operator=(ObjectiveFunction const& objective_function);

        Mtrx get_value(Mtrx const& x) const;
        void add_column(long double value);
        void log() const;

       private:
        ObjectiveFunction_T type;
        Mtrx                coefficients;
    };
    class LinearProblem {
       public:
        ~LinearProblem();
        LinearProblem(char const* _name);
        LinearProblem(const LinearProblem& linear_problem);
        LinearProblem& operator=(LinearProblem const& linear_problem);


        void add_variable(Variable* var);
        void load_problem(std::ifstream& file);
        void add_constraint(Constraint const& constraint);
        void set_objective_function(ObjectiveFunction const& objective_function);

        void log() const;
        void plot() const;
        void print_solution() const;
        void print_tableu_itr(Mtrx tableau, size_t itr) const;

        void solve();
        bool simplex();
        bool two_phase();
        void branch_and_bound();
        void process_to_standard_form();

       protected:
        ssize_t is_optimal(Mtrx const& tableau, ObjectiveFunction_T type, bool ignore,
                           ColumnSet const& _base);            // Checar esta funcion
        ssize_t min_ratio(Mtrx const& tableau, size_t index);  // Tambien revisar esta

        // Datos del problema lineal
        std::string             name;
        size_t                  solution_dimension;
        size_t                  original_solution_dimension;
        ObjectiveFunction       objective_function;
        ObjectiveFunction       aux_objective_function;
        ObjectiveFunction       plt_objctv_fnctn;  // Copia para poder graficar
        std::vector<Variable*>  variables;
        std::vector<Constraint> constraints;
        std::vector<Constraint> plt_cnstrnts;  // Copia para poder graficar
        std::vector<Constraint> binary_constraints;
        std::vector<Constraint> integer_constraints;
        std::vector<Constraint> no_negative_constraints;
        bool                    changed_sign         = false;
        bool                    binary_problem       = true;
        bool                    integer_problem      = false;
        bool                    artificial_constrait = false;
        // binary_problem (0-1 IP) es verdadero hasta que se demuestre lo contrario

        // Resultados
        Mtrx        solution;
        bool        feasible = false;
        long double solution_value;
    };

}  // namespace optimization

#endif
