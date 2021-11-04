#include "simplex.h"

namespace optimization {
    void ColumnSet::insert(size_t column) { columns.push_back(column); }

    void ColumnSet::remove(size_t column) {
        std::vector<size_t>::iterator itr = find(columns.begin(), columns.end(), column);
        if (itr != columns.end()) {
            columns.erase(itr);
        }
    }

    size_t& ColumnSet::operator[](size_t idx) { return columns[idx]; }

    void ColumnSet::log() const {
        for (std::vector<size_t>::const_iterator it = columns.begin(); it != columns.end(); it++) {
            std::cout << *it << " ";
        }
    }

    bool ColumnSet::contains(size_t column) const {
        return (find(columns.begin(), columns.end(), column) != columns.end());
    }

    size_t ColumnSet::size() const { return columns.size(); }
}  // namespace optimization
