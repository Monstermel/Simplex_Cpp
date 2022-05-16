#include <iostream>

#include "simplex.h"
namespace optimization {

    Constraint::Constraint(Mtrx const& cffcnts, Constraint_T _type, long double _value)
        : coefficients(cffcnts), type(_type), value(_value) {}

    int Constraint::size() const { return coefficients.cols(); }

    void Constraint::log() const {
        size_t aux = coefficients.cols();
        for (size_t i = 0; i < aux; i++) {
            std::cout << coefficients(i) << "\t";
        }

        switch (type) {
            case EQUAL:
                std::cout << "=\t";
                break;

            case LESS_EQUAL:
                std::cout << "<=\t";
                break;

            case GREATER_EQUAL:
                std::cout << ">=\t";
                break;

            default:
                break;
        }

        if (type == NOT_NEGATIVE || type == INTEGER || type == BINARY) {
            std::cout << std::endl;
        } else {
            std::cout << value << std::endl;
        }
    }

    void Constraint::add_column(long double _value) {
        size_t aux = coefficients.cols();
        Mtrx   row(1, aux + 1);
        for (size_t i = 0; i < aux; i++) {
            row(i) = coefficients(i);
        }
        row(aux)     = _value;
        coefficients = row;
    }

}  // namespace optimization
