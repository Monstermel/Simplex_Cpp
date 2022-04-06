
#include "simplex.h"

namespace optimization {

    ObjectiveFunction::ObjectiveFunction() {}

    ObjectiveFunction::ObjectiveFunction(ObjectiveFunction_T _type, Mtrx const& cffcnts)
        : type(_type), coefficients(cffcnts) {}

    ObjectiveFunction::ObjectiveFunction(const ObjectiveFunction& objective_function)
        : type(objective_function.type), coefficients(objective_function.coefficients) {}

    ObjectiveFunction& ObjectiveFunction::operator=(ObjectiveFunction const& objective_function) {
        type         = objective_function.type;
        coefficients = objective_function.coefficients;
        return *this;
    }

    void ObjectiveFunction::log() const {
        if (type == MIN) {
            std::cout << "MIN (";
        } else {
            std::cout << "MAX (";
        }

        size_t aux = coefficients.cols();
        for (size_t i = 0; i < aux; i++) {
            if (i != aux - 1) {
                std::cout << coefficients(i) << "  ";
            } else {
                std::cout << coefficients(i);
            }
        }
        std::cout << ") * X " << std::endl;
    }

    Mtrx ObjectiveFunction::get_value(Mtrx const& x) const { return coefficients * x; }

    void ObjectiveFunction::add_column(long double value) {
        size_t aux = coefficients.cols();
        Mtrx   row(1, aux + 1);
        for (size_t i = 0; i < aux; ++i) {
            row(i) = coefficients(i);
        }

        row(aux)     = value;
        coefficients = row;
    }

}  // namespace optimization
