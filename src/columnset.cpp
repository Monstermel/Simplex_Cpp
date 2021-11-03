#include "simplex.h"

namespace optimization {
    void ColumnSet::insert(size_t column) { columns.push_back(column); }

    void ColumnSet::substitute(size_t old_column, size_t new_column) {
        std::vector<size_t>::iterator itr = find(columns.begin(), columns.end(), old_column);
        if (itr != columns.end()) {
            *(itr) = new_column;
        }
    }

    void ColumnSet::remove(size_t column) {
        std::vector<size_t>::iterator itr = find(columns.begin(), columns.end(), column);
        if (itr != columns.end()) {
            columns.erase(itr);
        }
    }

    int ColumnSet::index_of(size_t column) const {
        size_t pos = 0;
        for (std::vector<size_t>::const_iterator it = columns.begin(); it != columns.end(); it++) {
            if (*(it) != column) {
                pos++;
            } else {
                return pos;
            }
        }
        return -1;
    }

    void ColumnSet::log() const {
        for (std::vector<size_t>::const_iterator it = columns.begin(); it != columns.end(); it++) {
            std::cout << *it << " ";
        }
    }

    bool ColumnSet::contains(size_t column) const {
        return (find(columns.begin(), columns.end(), column) != columns.end());
    }

    size_t& ColumnSet::column(size_t idx) { return columns.at(idx); }

    size_t ColumnSet::size() const { return columns.size(); }
}  // namespace optimization
