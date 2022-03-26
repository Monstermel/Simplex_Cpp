
#include "variable.h"

namespace optimization {

    Variable::Variable(LinealProblem* crtr, char const* _name) : creator(crtr), name(_name) {}
    Variable::~Variable() {}
    void Variable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx);
    }

    SplittedVariable::SplittedVariable(LinealProblem* crtr, char const* _name, AuxiliaryVariable* _aux)
        : Variable(crtr, _name), aux(_aux) {}
    SplittedVariable::~SplittedVariable() {}
    void SplittedVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx) - calculated_solution(aux->idx);
    }

    SlackVariable::SlackVariable(LinealProblem* crtr, char const* _name) : Variable(crtr, _name) {}
    SlackVariable::~SlackVariable() {}
    void SlackVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        (void) calculated_solution;
        (void) solution;
        (void) _idx;
    }

    AuxiliaryVariable::AuxiliaryVariable(LinealProblem* crtr, char const* _name, size_t _idx)
        : Variable(crtr, _name), idx(_idx) {}
    AuxiliaryVariable::~AuxiliaryVariable() {}
    void AuxiliaryVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        (void) calculated_solution;
        (void) solution;
        (void) _idx;
    }

}  // namespace optimization
