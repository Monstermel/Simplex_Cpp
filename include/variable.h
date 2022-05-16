#ifndef VARIABLE_H
#define VARIABLE_H

#include "simplex.h"

namespace optimization {

    class AuxiliaryVariable;

    class Variable {
        friend class LinearProblem;

       public:
        Variable(LinearProblem* crtr, char const* _name);
        virtual ~Variable();
        virtual void      process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
        virtual Variable* clone();

       protected:
        LinearProblem* creator;
        std::string    name;
        Variable_T     type;
    };

    class SplittedVariable : public Variable {
        friend class LinearProblem;

       public:
        SplittedVariable(LinearProblem* crtr, char const* _name, AuxiliaryVariable* _aux);
        ~SplittedVariable();
        void                      process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
        virtual SplittedVariable* clone();


       private:
        AuxiliaryVariable* aux;
    };

    class SlackVariable : public Variable {
        friend class LinearProblem;

       public:
        SlackVariable(LinearProblem* crtr, char const* _name);
        ~SlackVariable();
        void                   process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
        virtual SlackVariable* clone();
    };

    class ExcessVariable : public Variable {
        friend class LinearProblem;

       public:
        ExcessVariable(LinearProblem* crtr, char const* _name);
        ~ExcessVariable();
        void                    process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
        virtual ExcessVariable* clone();
    };

    class ArtificialVariable : public Variable {
        friend class LinearProblem;

       public:
        ArtificialVariable(LinearProblem* crtr, char const* _name);
        ~ArtificialVariable();
        void                        process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
        virtual ArtificialVariable* clone();
    };

    class AuxiliaryVariable : public Variable {
        friend class LinearProblem;
        friend class SplittedVariable;

       public:
        AuxiliaryVariable(LinearProblem* crtr, char const* _name, size_t _idx);
        ~AuxiliaryVariable();
        void                       process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx);
        virtual AuxiliaryVariable* clone();

       private:
        size_t idx;
    };
}  // namespace optimization

#endif
