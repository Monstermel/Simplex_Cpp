#ifndef VARIABLE_H
#define VARIABLE_H

#include "simplex.h"

namespace optimization {

    class AuxiliaryVariable;

    class Variable {
        friend class Simplex;

       public:
        Variable(Simplex* crtr, char const* _name);
        virtual ~Variable();
        virtual void process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);

       protected:
        Simplex*    creator;
        std::string name;
    };

    class SplittedVariable : public Variable {
        friend class Simplex;

       public:
        SplittedVariable(Simplex* crtr, char const* _name, AuxiliaryVariable* _aux);
        ~SplittedVariable();
        void process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);

       private:
        AuxiliaryVariable* aux;
    };

    class SlackVariable : public Variable {
        friend class Simplex;

       public:
        SlackVariable(Simplex* crtr, char const* _name);
        ~SlackVariable();
        void process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
    };

    class AuxiliaryVariable : public Variable {
        friend class Simplex;
        friend class SplittedVariable;

       public:
        AuxiliaryVariable(Simplex* crtr, char const* _name, size_t _idx);
        ~AuxiliaryVariable();
        void process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);

       private:
        size_t idx;
    };
}  // namespace optimization

#endif
