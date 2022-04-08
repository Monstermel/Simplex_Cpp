
#include "variable.h"

namespace optimization {

    Variable::Variable(LinearProblem* crtr, char const* _name) : creator(crtr), name(_name) {
        type = ORDINARY;
    }
    Variable::~Variable() {}
    void Variable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx);
    }
    Variable* Variable::clone() { return (new Variable(*this)); }  // Devuelve una copia de si mismo

    SplittedVariable::SplittedVariable(LinearProblem* crtr, char const* _name,
                                       AuxiliaryVariable* _aux)
        : Variable(crtr, _name), aux(_aux) {
        type = SPLITTED;
    }
    SplittedVariable::~SplittedVariable() {}
    void SplittedVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx) - calculated_solution(aux->idx);
    }
    SplittedVariable* SplittedVariable::clone() {
        return (new SplittedVariable(*this));
    }  // Devuelve una copia de si mismo

    SlackVariable::SlackVariable(LinearProblem* crtr, char const* _name) : Variable(crtr, _name) {
        type = SLACK;
    }
    SlackVariable::~SlackVariable() {}
    void SlackVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx);
    }
    SlackVariable* SlackVariable::clone() {
        return (new SlackVariable(*this));
    }  // Devuelve una copia de si mismo

    ExcessVariable::ExcessVariable(LinearProblem* crtr, char const* _name) : Variable(crtr, _name) {
        type = EXCESS;
    }
    ExcessVariable::~ExcessVariable() {}
    void ExcessVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx);
    }
    ExcessVariable* ExcessVariable::clone() {
        return (new ExcessVariable(*this));
    }  // Devuelve una copia de si mismo

    ArtificialVariable::ArtificialVariable(LinearProblem* crtr, char const* _name)
        : Variable(crtr, _name) {
        type = ARTIFICIAL;
    }
    ArtificialVariable::~ArtificialVariable() {}
    void ArtificialVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        solution(_idx) = calculated_solution(_idx);
    }
    ArtificialVariable* ArtificialVariable::clone() {
        return (new ArtificialVariable(*this));
    }  // Devuelve una copia de si mismo

    AuxiliaryVariable::AuxiliaryVariable(LinearProblem* crtr, char const* _name, size_t _idx)
        : Variable(crtr, _name), idx(_idx) {
        type = AUXILIARY;
    }
    AuxiliaryVariable::~AuxiliaryVariable() {}
    void AuxiliaryVariable::process(Mtrx& calculated_solution, Mtrx& solution, size_t _idx) {
        (void) calculated_solution;
        (void) solution;
        (void) _idx;
    }
    AuxiliaryVariable* AuxiliaryVariable::clone() {
        return (new AuxiliaryVariable(*this));
    }  // Devuelve una copia de si mismo

}  // namespace optimization
