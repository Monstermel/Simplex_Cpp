#include "graph.h"

#include <fstream>
#include <iostream>
#include <queue>

#include "config.h"
#include "simplex.h"


namespace optimization {
    GraphProblem::GraphProblem(char const* _name) : name(_name) {}

    void GraphProblem::load_problem(std::ifstream& file) {
        // Variable enum para determiar que bloque estamos leyendo
        ParsingContext current_parsing_block = NONE;
        // Variable que cuenta el numero de nodos registradas
        size_t current_node = 0;

        while (!file.eof()) {
            // Extraemos una linea hasta encontrar \n
            std::string buffer_init;
            getline(file, buffer_init);
            std::stringstream buffer(buffer_init);
            std::string       token;
            buffer >> token;
            if (token.length()) {
                if (token == "Nodos")
                    current_parsing_block = VARS;
                else if (token == "Aristas")
                    current_parsing_block = CONSTRAINTS;
                else if (token == "Objetivo")
                    current_parsing_block = OBJECTIVE;
                else {
                    if (current_parsing_block == NONE) {
                        throw("Indentificador de bloque invalido: Ln 32:graph.cpp");
                    }
                    switch (current_parsing_block) {
                        case VARS: {
                            std::string node_name;
                            buffer >> node_name;
                            std::string value;
                            buffer >> value;
                            if (value.empty()) {
                                value = "0";
                            }
                            long long b_current_node;
                            try {
                                b_current_node = std::stoll(value);
                            } catch (...) {
                                throw("Definicion de nodos invalida: Ln 47:graph.cpp");
                            }
                            if (b_current_node < 0) {
                                nodes.emplace(node_name, current_node, DEMAND);
                                num_demand_pnts++;
                            } else if (b_current_node > 0) {
                                nodes.emplace(node_name, current_node, SUPPLY);
                                num_supply_pnts++;
                            } else {
                                nodes.emplace(node_name, current_node, DEFAULT);
                                num_supply_pnts++;
                                num_demand_pnts++;
                            }
                            b_values.push_back(b_current_node);
                            current_node++;
                        } break;

                        case CONSTRAINTS: {
                            auto itr = nodes.find(token);
                            if (itr == nodes.end()) {
                                throw("Definicion de aristas invalida: Ln 67:graph.cpp");
                            }
                            if (adjacency.cols() == 0) {
                                adjacency = MtrxI::Zero(nodes.size(), nodes.size());
                                adjacency.fill(INF);
                            }
                            // Por ahora solo voy a dejar el token > para marcar aristas
                            // dirigidas
                            buffer >> token;
                            // Procesado nodo:peso_arista
                            // getline no quita espacios en blanco
                            while (!buffer.eof()) {
                                std::getline(buffer, token, ':');
                                token.erase(remove_if(token.begin(), token.end(), isspace),
                                            token.end());
                                if (token.empty()) {
                                    break;
                                }
                                auto itr_1 = nodes.find(token);
                                if (itr_1 == nodes.end()) {
                                    throw("Definicion de aristas invalida: Ln 87:graph.cpp");
                                }
                                buffer >> token;
                                try {
                                    adjacency(itr->index, itr_1->index) = std::stoll(token);
                                } catch (...) {
                                    throw("Definicion de aristas invalida: Ln 93:graph.cpp");
                                }
                            }
                        } break;
                        case OBJECTIVE: {
                            if (token == "TRANSPORTATION") {
                                type = TRANSPORTATION;
                            } else if (token == "ASSIGNMENT") {
                                type = ASSIGNMENT;
                            } else if (token == "MIN_SPANNING_TREE") {
                                type = MIN_SPANNING_TREE;
                            } else if (token == "SHORTEST_PATH") {
                                type = SHORTEST_PATH_FLOYD;
                                while (!buffer.eof()) {
                                    std::string node_name;
                                    buffer >> node_name;
                                    auto itr = nodes.find(node_name);
                                    if (itr == nodes.end()) {
                                        throw(
                                            "Definicion de objetivo invalida: Ln "
                                            "129:graph.cpp");
                                    }
                                    source = itr->index;
                                    type   = SHORTEST_PATH_DIJKSTRA;
                                }

                            } else if (token == "MAX_FLOW") {
                                type = MAX_FLOW;
                                std::string node_name;
                                // SOURCE
                                buffer >> node_name;
                                auto itr = nodes.find(node_name);
                                if (itr == nodes.end()) {
                                    throw("Definicion de objetivo invalida: Ln 142:graph.cpp");
                                }
                                source = itr->index;

                                buffer >> node_name;  // descartamos >

                                // DESTINATION
                                buffer >> node_name;
                                itr = nodes.find(node_name);
                                if (itr == nodes.end()) {
                                    throw("Definicion de objetivo invalida: Ln 152:graph.cpp");
                                }
                                destination = itr->index;

                            } else {
                                throw("Definicion de objetivo invalida Ln 157:graph.cpp");
                            }
                        }
                        default:
                            break;
                    }
                }
            }
        }
        file.close();
        assert(type != NOT_ASSIGNED);

        return;
    }

    std::string GraphProblem::get_node_name(ssize_t idx) const {
        for (auto nodo : nodes) {
            if (nodo.index == idx) {
                return nodo.name;
            }
        }
        return {};
    }

    int GraphProblem::is_not_balanced() {
        total_demand = 0;
        total_supply = 0;
        for (size_t i = 0; i < b_values.size(); i++) {
            if (b_values[i] < 0) {
                total_demand += b_values[i];
            } else if (b_values[i] > 0) {
                total_supply += b_values[i];
            }
        }
        if (total_supply == -(total_demand))
            return 0;
        else if (total_supply < -(total_demand))
            return -1;
        else
            return 1;
    }

    void GraphProblem::add_dummy_node() {
        adjacency.conservativeResize(nodes.size() + 1, nodes.size() + 1);
        adjacency(nodes.size(), nodes.size()) = INF;
        for (size_t i = 0; i < nodes.size(); i++) {
            adjacency(i, nodes.size()) = 0;
            adjacency(nodes.size(), i) = INF;
        }
        b_values.push_back(-(total_supply + total_demand));
        nodes.emplace("dummy", nodes.size(), DEMAND);
        num_demand_pnts++;
        // Ahora esta balanceado
        total_demand = -total_supply;
        return;
    }

    void GraphProblem::build_tableau(MtrxI& costs, MtrxI& tableau,
                                     std::vector<ssize_t>& supply_relation,
                                     std::vector<ssize_t>& demand_relation) const {
        costs   = MtrxI::Zero(num_supply_pnts, num_demand_pnts);
        tableau = MtrxI::Zero(num_supply_pnts + 1, num_demand_pnts + 1);
        // Tomamos cada nodo y lo relacionamos con una fila o columna del tanleau
        for (auto node : nodes) {
            if (node.type == SUPPLY) {
                supply_relation.push_back(node.index);
            } else if (node.type == DEMAND) {
                demand_relation.push_back(node.index);
            } else if (node.type == DEFAULT) {
                supply_relation.push_back(node.index);
                demand_relation.push_back(node.index);
            }
        }
        assert(supply_relation.size() == num_supply_pnts &&
               demand_relation.size() == num_demand_pnts);
        assert(total_supply == -total_demand);
        // Cargamos la supply y demand: faltan los de transporte: ERROR
        for (ssize_t i = 0; i < tableau.rows() - 1; i++) {
            if (!(b_values[supply_relation[i]])) {
                tableau(i, tableau.cols() - 1) = total_supply;
            } else {
                tableau(i, tableau.cols() - 1) = b_values[supply_relation[i]];
            }
        }
        for (ssize_t i = 0; i < tableau.cols() - 1; i++) {
            if (!(b_values[demand_relation[i]])) {
                tableau(tableau.rows() - 1, i) = total_supply;
            } else {
                tableau(tableau.rows() - 1, i) = -b_values[demand_relation[i]];
            }
        }
        // Cargamos los costos
        for (ssize_t i = 0; i < costs.rows(); i++) {
            for (ssize_t j = 0; j < costs.cols(); j++) {
                costs(i, j) = adjacency(supply_relation[i], demand_relation[j]);
            }
        }
        return;
    }

    bool GraphProblem::bfs_not_found(std::vector<bool> const& supply_available,
                                     std::vector<bool> const& demand_available,
                                     size_t&                  total_supply_available,
                                     size_t&                  total_demand_available) const {
        total_supply_available = 0;
        for (bool foo : supply_available) {
            if (foo) {
                total_supply_available++;
            }
        }
        total_demand_available = 0;
        for (bool foo : demand_available) {
            if (foo) {
                total_demand_available++;
            }
        }
        if (total_demand_available + total_supply_available > 1) {
            return true;
        }
        return false;
    }

    void GraphProblem::vogel_method(MtrxI& costs, MtrxI& tableau, MtrxB& basic_variables) const {
        std::vector<bool> supply_available(num_supply_pnts, true);
        std::vector<bool> demand_available(num_demand_pnts, true);
        size_t            total_supply_available;
        size_t            total_demand_available;
        ssize_t           s_num_supply_pnts = num_supply_pnts;
        ssize_t           s_num_demand_pnts = num_demand_pnts;
        while (bfs_not_found(supply_available, demand_available, total_supply_available,
                             total_demand_available)) {
            // Buscamos los dos costos mas pequeños por columna y fila
            std::vector<pair_idx> min_idx_row;
            for (ssize_t i = 0; i < s_num_supply_pnts; i++) {
                min_idx_row.push_back({-1, -1});
                if (!supply_available[i]) {
                    continue;
                }
                bool found = false;
                for (ssize_t j = 0; j < s_num_demand_pnts; j++) {
                    if (!demand_available[j]) {
                        continue;
                    }
                    if (!found) {
                        min_idx_row[i].first = j;
                        found                = true;
                        if (total_demand_available == 1) {
                            min_idx_row[i].second = j;
                            break;
                        }
                    } else if (costs(i, j) <= costs(i, min_idx_row[i].first)) {
                        min_idx_row[i].first = j;
                    }
                }
                found = false;
                for (ssize_t j = 0; j < s_num_demand_pnts; j++) {
                    if (total_demand_available == 1) {
                        break;
                    }
                    if (!demand_available[j]) {
                        continue;
                    }
                    if (!found) {
                        min_idx_row[i].second = j;
                        found                 = true;
                    } else if (costs(i, j) <= costs(i, min_idx_row[i].second) &&
                               costs(i, min_idx_row[i].first) <= costs(i, j) &&
                               min_idx_row[i].first != j) {
                        min_idx_row[i].second = j;
                    }
                }
            }
            std::vector<pair_idx> min_idx_col;
            for (ssize_t i = 0; i < s_num_demand_pnts; i++) {
                min_idx_col.push_back({-1, -1});
                if (!demand_available[i]) {
                    continue;
                }
                bool found = false;
                for (ssize_t j = 0; j < s_num_supply_pnts; j++) {
                    if (!supply_available[j]) {
                        continue;
                    }
                    if (!found) {
                        min_idx_col[i].first = j;
                        found                = true;
                        if (total_supply_available == 1) {
                            min_idx_col[i].second = j;
                            break;
                        }
                    } else if (costs(j, i) <= costs(min_idx_col[i].first, i)) {
                        min_idx_col[i].first = j;
                    }
                }
                found = false;
                for (ssize_t j = 0; j < s_num_supply_pnts; j++) {
                    if (total_supply_available == 1) {
                        break;
                    }
                    if (!supply_available[j]) {
                        continue;
                    }
                    if (!found) {
                        min_idx_col[i].second = j;
                        found                 = true;
                    } else if (costs(j, i) <= costs(min_idx_col[i].second, i) &&
                               costs(min_idx_col[i].first, i) <= costs(j, i) &&
                               min_idx_col[i].first != j) {
                        min_idx_col[i].second = j;
                    }
                }
            }
            // Buscamos la penalizacion maxima
            long long max_penalty{0};
            pair_idx  penalty_idx;
            for (ssize_t i = 0; i < costs.rows(); i++) {
                if (supply_available[i]) {
                    max_penalty =
                        std::abs(costs(i, min_idx_row[i].first) - costs(i, min_idx_row[i].second));
                    penalty_idx = {i, min_idx_row[i].first};
                    break;
                }
            }
            // Por filas
            for (ssize_t i = 0; i < costs.rows(); i++) {
                if (!supply_available[i]) {
                    continue;
                }
                if (!tableau(i, tableau.cols() - 1)) {
                    max_penalty = std::numeric_limits<long long>::max();
                    penalty_idx = {i, min_idx_row[i].first};
                } else if (max_penalty < std::abs(costs(i, min_idx_row[i].first) -
                                                  costs(i, min_idx_row[i].second))) {
                    max_penalty =
                        std::abs(costs(i, min_idx_row[i].first) - costs(i, min_idx_row[i].second));
                    penalty_idx = {i, min_idx_row[i].first};
                }
            }
            // Por columnas
            for (ssize_t i = 0; i < costs.cols(); i++) {
                if (!demand_available[i]) {
                    continue;
                }
                if (!tableau(tableau.rows() - 1, i)) {
                    max_penalty = std::numeric_limits<long long>::max();
                    penalty_idx = {min_idx_col[i].first, i};
                } else if (max_penalty < std::abs(costs(min_idx_col[i].first, i) -
                                                  costs(min_idx_col[i].second, i))) {
                    max_penalty =
                        std::abs(costs(min_idx_col[i].first, i) - costs(min_idx_col[i].second, i));
                    penalty_idx = {min_idx_col[i].first, i};
                }
            }
            // Elegimos la variable basica
            basic_variables(penalty_idx.first, penalty_idx.second) = true;
            tableau(penalty_idx.first, penalty_idx.second) =
                std::min(tableau(penalty_idx.first, tableau.cols() - 1),
                         tableau(tableau.rows() - 1, penalty_idx.second));
            tableau(penalty_idx.first, tableau.cols() - 1) -=
                tableau(penalty_idx.first, penalty_idx.second);
            tableau(tableau.rows() - 1, penalty_idx.second) -=
                tableau(penalty_idx.first, penalty_idx.second);
            if (!tableau(penalty_idx.first, tableau.cols() - 1)) {
                supply_available[penalty_idx.first] = false;
            } else if (!tableau(tableau.rows() - 1, penalty_idx.second)) {
                demand_available[penalty_idx.second] = false;
            }
        }
    }

    pair_idx GraphProblem::is_optimal(MtrxI const& costs, MtrxB const& basic_variables) const {
        const long long        NOT_SET = std::numeric_limits<long long>::max();
        std::vector<long long> U_values(num_supply_pnts, NOT_SET);
        std::vector<long long> V_values(num_demand_pnts, NOT_SET);
        std::queue<pair_idx>   not_solved_yet;
        U_values[0] = 0;
        for (ssize_t i = 0; i < basic_variables.rows(); i++) {
            for (ssize_t j = 0; j < basic_variables.cols(); j++) {
                if (basic_variables(i, j)) {
                    if (U_values[i] == NOT_SET && V_values[j] == NOT_SET) {
                        not_solved_yet.push({i, j});
                    } else if (U_values[i] == NOT_SET) {
                        U_values[i] = costs(i, j) - V_values[j];
                    } else if (V_values[j] == NOT_SET) {
                        V_values[j] = costs(i, j) - U_values[i];
                    }
                }
            }
        }
        while (!not_solved_yet.empty()) {
            pair_idx current = not_solved_yet.front();
            not_solved_yet.pop();
            if (U_values[current.first] == NOT_SET && V_values[current.second] == NOT_SET) {
                not_solved_yet.push({current.first, current.second});
            } else if (U_values[current.first] == NOT_SET) {
                U_values[current.first] =
                    costs(current.first, current.second) - V_values[current.second];
            } else if (V_values[current.second] == NOT_SET) {
                V_values[current.second] =
                    costs(current.first, current.second) - U_values[current.first];
            }
        }
        MtrxI    shadow_prices = MtrxI::Zero(num_supply_pnts, num_demand_pnts);
        pair_idx optimal_test  = {0, 0};
        // Calculamos los precios sombra
        for (ssize_t i = 0; i < shadow_prices.rows(); i++) {
            for (ssize_t j = 0; j < shadow_prices.cols(); j++) {
                shadow_prices(i, j) = U_values[i] + V_values[j] - costs(i, j);
                optimal_test =
                    (shadow_prices(i, j) < shadow_prices(optimal_test.first, optimal_test.second))
                        ? optimal_test
                        : pair_idx{i, j};
            }
        }
        if (0 < shadow_prices(optimal_test.first, optimal_test.second)) {
            return {optimal_test.first, optimal_test.second};
        }
        return {-1, -1};
    }

    void GraphProblem::find_loop(pair_idx const& idx, MtrxB const& basic_variables,
                                 std::vector<pair_idx>& loop) const {
        assert(idx.first != -1 && idx.second != -1);
        typedef std::pair<pair_idx, pair_idx> Pair_Pairs;
        typedef std::vector<Pair_Pairs>       Vector_Pairs;
        typedef std::vector<Vector_Pairs>     Vector_Vector_Pairs;
        MtrxI visited                  = MtrxI::Zero(num_supply_pnts, num_demand_pnts);
        visited(idx.first, idx.second) = 1;
        std::queue<pair_idx> Queue;
        Queue.push(idx);
        Vector_Vector_Pairs parents;
        parents.push_back({{idx, {-1, -1}}});
        pair_idx end_cycle;
        while (!Queue.empty()) {
            pair_idx current = Queue.front();
            Queue.pop();
            if (static_cast<long long>(parents.size()) !=
                visited(current.first, current.second) + 1) {
                parents.push_back(Vector_Pairs());
            }
            for (ssize_t i = 0; i < basic_variables.rows(); i++) {
                if (i == current.first) {
                    continue;
                }
                if (visited(i, current.second) &&
                    visited(i, current.second) == (visited(current.first, current.second) + 1)) {
                    parents[visited(i, current.second) - 1].push_back(
                        {{i, current.second}, current});
                    end_cycle = {i, current.second};
                    goto found;
                }
                if (basic_variables(i, current.second) && !visited(i, current.second)) {
                    visited(i, current.second) = visited(current.first, current.second) + 1;
                    parents[visited(i, current.second) - 1].push_back(
                        {{i, current.second}, current});
                    Queue.push({i, current.second});
                }
            }
            for (ssize_t i = 0; i < basic_variables.cols(); i++) {
                if (i == current.second) {
                    continue;
                }
                if (visited(current.first, i) &&
                    visited(current.first, i) == (visited(current.first, current.second) + 1)) {
                    parents[visited(current.first, i) - 1].push_back({{current.first, i}, current});
                    end_cycle = {current.first, i};
                    goto found;
                }
                if (basic_variables(current.first, i) && !visited(current.first, i)) {
                    visited(current.first, i) = visited(current.first, current.second) + 1;
                    parents[visited(current.first, i) - 1].push_back({{current.first, i}, current});
                    Queue.push({current.first, i});
                }
            }
        }
    found:
        // Reconstruimos el loop
        std::vector<pair_idx> first_half;
        std::vector<pair_idx> second_half;
        for (auto a : parents[parents.size() - 1]) {
            if (a.first.first == end_cycle.first && a.first.second == end_cycle.second &&
                first_half.empty()) {
                first_half.push_back(a.first);
                first_half.push_back(a.second);
            } else if (a.first.first == end_cycle.first && a.first.second == end_cycle.second) {
                second_half.push_back(a.second);
            }
            if (!first_half.empty() && !second_half.empty()) {
                break;
            }
        }
        for (ssize_t i = parents.size() - 2; i; i--) {
            for (auto obj : parents[i]) {
                if (obj.first.first == first_half.back().first &&
                    obj.first.second == first_half.back().second) {
                    first_half.push_back(obj.second);
                } else if (obj.first.first == second_half.back().first &&
                           obj.first.second == second_half.back().second) {
                    second_half.push_back(obj.second);
                }
            }
        }
        first_half.pop_back();
        loop.insert(loop.end(), second_half.rbegin(), second_half.rend());
        loop.insert(loop.end(), first_half.begin(), first_half.end());
        return;
    }

    void GraphProblem::transportation() {
        switch (is_not_balanced()) {
            case -1:
                throw("Problema infactible");
                break;
            case 1:
                add_dummy_node();
                break;
            default:
                break;
        }
        MtrxI                costs;
        MtrxI                tableau;
        std::vector<ssize_t> supply_relation;
        std::vector<ssize_t> demand_relation;
        build_tableau(costs, tableau, supply_relation, demand_relation);
        MtrxB basic_variables = MtrxB::Zero(num_supply_pnts, num_demand_pnts);
        vogel_method(costs, tableau, basic_variables);
        for (pair_idx enter_var = is_optimal(costs, basic_variables); enter_var.first != -1;
             enter_var          = is_optimal(costs, basic_variables)) {
            std::vector<pair_idx> loop;
            find_loop(enter_var, basic_variables, loop);
            pair_idx min_obv = loop[1];
            for (size_t i = 3; i < loop.size(); i += 2) {
                min_obv = (tableau(min_obv.first, min_obv.second) <
                           tableau(loop[i].first, loop[i].second))
                              ? min_obv
                              : loop[i];
            }
            const long long min_obs_value = tableau(min_obv.first, min_obv.second);
            for (size_t i = 0; i < loop.size(); i++) {
                if (!(i % 2)) {
                    tableau(loop[i].first, loop[i].second) += min_obs_value;
                } else {
                    tableau(loop[i].first, loop[i].second) -= min_obs_value;
                }
            }
            basic_variables(min_obv.first, min_obv.second)     = false;
            basic_variables(enter_var.first, enter_var.second) = true;
        }
        std::cout << "Fo = " << tableau(tableau.rows() - 1, tableau.cols() - 1) << '\n';
        for (ssize_t i = 0; i < tableau.rows() - 1; i++) {
            for (ssize_t j = 0; j < tableau.cols() - 1; j++) {
                // antes estaba tableau(i, j), verificar que esta mal eso
                if (basic_variables(i, j)) {
                    std::cout << get_node_name(supply_relation[i]);
                    std::cout << "->";
                    std::cout << get_node_name(demand_relation[j]);
                    std::cout << ": " << tableau(i, j) << '\n';
                }
            }
        }
    }

    void GraphProblem::build_hungarian_tableau(MtrxI& costs, std::vector<ssize_t>& supply_relation,
                                               std::vector<ssize_t>& demand_relation) const {
        // Tomamos cada nodo y lo relacionamos con una fila o columna del tanleau
        for (auto node : nodes) {
            if (node.type == SUPPLY) {
                supply_relation.push_back(node.index);
            } else if (node.type == DEMAND) {
                demand_relation.push_back(node.index);
            }
        }
        assert(supply_relation.size() == num_supply_pnts &&
               demand_relation.size() == num_demand_pnts);
        assert(total_supply == -total_demand);
        // Cargamos los costos
        for (ssize_t i = 0; i < costs.rows(); i++) {
            for (ssize_t j = 0; j < costs.cols(); j++) {
                costs(i, j) = adjacency(supply_relation[i], demand_relation[j]);
            }
        }
        return;
    }

    void GraphProblem::subtract_min(MtrxI& costs) const {
        for (ssize_t i = 0; i < costs.rows(); i++) {
            long long min = costs(i, 0);
            for (ssize_t j = 1; j < costs.cols(); j++) {
                if (costs(i, j) < min) {
                    min = costs(i, j);
                }
            }
            for (ssize_t j = 0; j < costs.cols(); j++) {
                costs(i, j) -= min;
            }
        }
        for (ssize_t j = 0; j < costs.cols(); j++) {
            long long min = costs(0, j);
            for (ssize_t i = 1; i < costs.rows(); i++) {
                if (costs(i, j) < min) {
                    min = costs(i, j);
                }
            }
            for (ssize_t i = 0; i < costs.rows(); i++) {
                costs(i, j) -= min;
            }
        }
    }

    bool GraphProblem::try_kuhn(ssize_t i, std::vector<bool>& used, std::vector<ssize_t>& match_col,
                                MtrxI const& costs) const {
        if (used[i]) {
            return false;
        }
        used[i] = true;
        for (ssize_t j = 0; j < costs.cols(); j++) {
            if (!costs(i, j)) {
                if (match_col[j] == -1 || try_kuhn(match_col[j], used, match_col, costs)) {
                    match_col[j] = i;
                    return true;
                }
            }
        }
        return false;
    }

    void GraphProblem::draw_lines(MtrxI const& costs, MtrxI& lines,
                                  std::vector<ssize_t> const& match_col) const {
        MtrxI zero_used{MtrxI::Constant(costs.rows(), costs.cols(), -1)};
        for (ssize_t i = 0; i < costs.rows(); i++) {
            for (ssize_t j = 0; j < costs.cols(); j++) {
                if (!costs(i, j)) {
                    zero_used(i, j) = 0;
                }
            }
        }

        for (ssize_t j = 0; j < costs.cols(); j++) {
            if (match_col[j] != -1) {
                size_t zeros_row = 0;
                for (ssize_t k = 0; k < zero_used.cols(); k++) {
                    if (!zero_used(match_col[j], k)) {
                        zeros_row++;
                    }
                }
                size_t zeros_col = 0;
                for (ssize_t k = 0; k < zero_used.rows(); k++) {
                    if (!zero_used(k, j)) {
                        zeros_col++;
                    }
                }
                if (zeros_col <= zeros_row) {
                    for (ssize_t k = 0; k < zero_used.cols(); k++) {
                        if (zero_used(match_col[j], k) != -1) {
                            zero_used(match_col[j], k)++;
                        }
                        lines(match_col[j], k)++;
                    }
                } else {
                    for (ssize_t k = 0; k < zero_used.rows(); k++) {
                        if (zero_used(k, j) != -1) {
                            zero_used(k, j)++;
                        }
                        lines(k, j)++;
                    }
                }
            }
        }
    }

    ssize_t GraphProblem::find_min_lines(MtrxI const& costs, MtrxI& lines) const {
        // Kuhn's Algorithm for Maximum Bipartite Matching
        std::vector<ssize_t> match_col;
        std::vector<bool>    used;
        match_col.assign(costs.rows(), -1);

        ssize_t num_lines = 0;
        for (ssize_t i = 0; i < costs.rows(); i++) {
            used.assign(costs.rows(), false);
            if (try_kuhn(i, used, match_col, costs)) {
                num_lines++;
            }
        }

        lines = MtrxI::Zero(costs.rows(), costs.cols());
        draw_lines(costs, lines, match_col);

        return num_lines;
    }

    /*
    Solucion correcta: 8 lineas
    costs = MtrxI{{0, 194, 112, 191, 157, 0, 115, 112, 199},
                  {113, 119, 117, 132, 142, 113, 0, 135, 116},
                  {109, 111, 131, 156, 138, 129, 116, 131, 0},
                  {181, 151, 139, 0, 110, 137, 124, 167, 140},
                  {194, 0, 134, 159, 123, 142, 127, 130, 111},
                  {171, 137, 139, 0, 0, 147, 132, 171, 148},
                  {171, 141, 143, 114, 0, 143, 128, 171, 144},
                  {180, 110, 0, 153, 137, 0, 113, 0, 197},
                  {0, 194, 0, 189, 157, 118, 121, 0, 105}};

    */

    void GraphProblem::hungarian_method() {
        switch (is_not_balanced()) {
            case -1:
                throw("Problema infactible");
                break;
            case 1:
                add_dummy_node();
                break;
            default:
                break;
        }
        std::vector<ssize_t> supply_relation;
        std::vector<ssize_t> demand_relation;
        MtrxI                costs(num_supply_pnts, num_demand_pnts);
        build_hungarian_tableau(costs, supply_relation, demand_relation);
        if (VERBOSE) {
            std::cout << "Costos iniciales:" << '\n'
                      << costs << "\n\n"
                      << "###################################################################\n";
        }

        MtrxI lines;
        MtrxI rcosts{costs};
        subtract_min(rcosts);
        for (ssize_t min_lines = find_min_lines(rcosts, lines); min_lines != rcosts.rows();
             min_lines         = find_min_lines(rcosts, lines)) {
            if (VERBOSE) {
                std::cout
                    << "Costos reducidos:\n"
                    << rcosts << "\n\n"
                    << "Lineas:\n"
                    << lines << "\n\n"
                    << "###################################################################\n";
            }

            long long min_element = -1;
            for (ssize_t i = 0; i < rcosts.rows(); i++) {
                for (ssize_t j = 0; j < rcosts.cols(); j++) {
                    if (min_element == -1 && !lines(i, j)) {
                        min_element = rcosts(i, j);
                    } else if (!lines(i, j) && rcosts(i, j) < min_element) {
                        min_element = rcosts(i, j);
                    }
                }
            }


            for (ssize_t i = 0; i < rcosts.rows(); i++) {
                for (ssize_t j = 0; j < rcosts.cols(); j++) {
                    if (!lines(i, j)) {
                        rcosts(i, j) -= min_element;
                    } else if (lines(i, j) == 2) {
                        rcosts(i, j) += min_element;
                    }
                }
            }
        }
        if (VERBOSE) {
            std::cout << "Costos reducidos:\n"
                      << rcosts << "\n\n"
                      << "Lineas:\n"
                      << lines << "\n\n"
                      << "###################################################################\n";
        }

        long long         FO       = 0;
        MtrxB             solution = MtrxB::Zero(costs.rows(), costs.cols());
        std::vector<bool> row_avalible(costs.rows(), true);
        std::vector<bool> col_avalible(costs.cols(), true);
        for (ssize_t i = 0; i < solution.rows(); i++) {
            for (ssize_t j = 0; j < solution.cols(); j++) {
                if (!rcosts(i, j) && row_avalible[i] && col_avalible[j]) {
                    FO += costs(i, j);
                    solution(i, j)  = true;
                    row_avalible[i] = false;
                    col_avalible[j] = false;
                }
            }
        }

        if (VERBOSE) {
            for (size_t i = 0; i < supply_relation.size(); i++) {
                std::cout << '[' << i << "][] = " << get_node_name(supply_relation[i]) << '\n';
            }
            for (size_t j = 0; j < demand_relation.size(); j++) {
                std::cout << "[][" << j << "] = " << get_node_name(demand_relation[j]) << '\n';
            }
            std::cout << '\n';
        }

        // RESULTADOS
        std::cout << "FO = " << FO << "\n";
        for (ssize_t i = 0; i < solution.rows(); i++) {
            for (ssize_t j = 0; j < solution.cols(); j++) {
                if (solution(i, j)) {
                    std::cout << get_node_name(supply_relation[i]) << "->"
                              << get_node_name(demand_relation[j]) << '\n';
                }
            }
        }
        ////////////
    }

    size_t GraphProblem::find_set(ssize_t v) {
        if (v == parent[v]) {
            return v;
        }
        return parent[v] = find_set(parent[v]);
    }

    void GraphProblem::union_sets(ssize_t a, ssize_t b) {
        a = find_set(a);
        b = find_set(b);
        if (a != b) {
            if (rank[a] < rank[b]) {
                std::swap(a, b);
            }
            parent[b] = a;
            if (rank[a] == rank[b]) {
                rank[a]++;
            }
        }
    }

    void GraphProblem::min_spanning_tree() {
        std::vector<Edge> edges;
        for (ssize_t i = 0; i < adjacency.rows(); i++) {
            for (ssize_t j = 0; j < adjacency.cols(); j++) {
                if (adjacency(i, j) != INF) {
                    edges.push_back({i, j, adjacency(i, j)});
                }
            }
        }

        std::vector<Edge> result;
        parent.resize(edges.size());
        rank.resize(edges.size());
        ssize_t edges_size = static_cast<ssize_t>(edges.size());
        for (ssize_t i = 0; i < edges_size; i++) {
            parent[i] = i;
            rank[i]   = 0;
        }

        std::sort(edges.begin(), edges.end());
        if (VERBOSE) {
            std::cout << "Aristas ordenadas:\n";
            for (auto edge : edges) {
                std::cout << '\t' << get_node_name(edge.u) << "->" << get_node_name(edge.v) << ": "
                          << edge.weight << '\n';
            }
            std::cout << "Verificacion anti ciclos:\n";
        }

        ssize_t cost = 0;
        for (Edge e : edges) {
            if (find_set(e.u) != find_set(e.v)) {
                if (VERBOSE) {
                    std::cout << '\t' << get_node_name(e.u) << " y " << get_node_name(e.v)
                              << " se pueden unir\n";
                }
                cost += e.weight;
                result.push_back(e);
                union_sets(e.u, e.v);
                if (static_cast<ssize_t>(result.size() + 1) == adjacency.rows()) {
                    break;
                }
            }
        }
        // RESULTADOS
        std::cout << "Costo: " << cost << '\n';
        for (Edge edge : result) {
            std::cout << get_node_name(edge.u) << "->" << get_node_name(edge.v) << ": "
                      << edge.weight << '\n';
        }
        /////////////
    }

    std::vector<ssize_t> GraphProblem::restore_path(ssize_t s, ssize_t t,
                                                    std::vector<ssize_t> const& p) const {
        std::vector<ssize_t> path;

        for (ssize_t v = t; v != s; v = p[v]) {
            path.push_back(v);
        }
        path.push_back(s);

        std::reverse(path.begin(), path.end());
        return path;
    }

    void GraphProblem::dijkstra_algorithm() {
        using std::vector, std::pair;
        vector<vector<pair<ssize_t, ssize_t>>> adj;
        for (ssize_t i = 0; i < adjacency.rows(); i++) {
            adj.push_back({});
            for (ssize_t j = 0; j < adjacency.cols(); j++) {
                if (adjacency(i, j) != INF) {
                    adj[i].push_back({j, adjacency(i, j)});
                }
            }
        }
        ssize_t         n = adj.size();
        vector<bool>    used(n, false);
        vector<ssize_t> distance(n, INF);
        vector<ssize_t> predecessors(n, -1);
        distance[source] = 0;

        for (ssize_t i = 0; i < n; i++) {
            ssize_t v = -1;
            for (ssize_t j = 0; j < n; j++) {
                if (!used[j] && (v == -1 || distance[j] < distance[v])) {
                    v = j;
                }
            }

            if (distance[v] == INF) {
                break;
            }

            std::string curent_name;

            if (VERBOSE) {
                curent_name = get_node_name(v);
                std::cout << "Nodo visitado: " << curent_name << '\n';
                std::cout << "Distancia minima: " << distance[v] << '\n';
                std::cout << "Nodos vecinos:\n";
            }

            used[v] = true;
            for (auto edge : adj[v]) {
                ssize_t to  = edge.first;
                ssize_t len = edge.second;
                if (VERBOSE) {
                    std::cout << "\tNodo vecino: " << get_node_name(to) << '\n';
                    std::cout << "\tDistancia desde " << curent_name << ": " << distance[v] + len
                              << '\n';
                }
                if (distance[v] + len < distance[to]) {
                    distance[to]     = distance[v] + len;
                    predecessors[to] = v;
                }
                if (VERBOSE) {
                    std::cout << "\tDistancia minima encontrada: " << distance[to] << "\n\n";
                }
            }
            if (VERBOSE) {
                std::cout
                    << "###################################################################\n";
            }
        }

        // RESULTADOS
        for (ssize_t i = 0; i < adjacency.rows(); i++) {
            if (i != source) {
                std::cout << distance[i] << ": ";
                for (auto itr : restore_path(source, i, predecessors)) {
                    std::cout << get_node_name(itr) << " ";
                }
                std::cout << "\n";
            }
        }
        /////////////
    }

    void GraphProblem::floyd_warshall_algorithm() {
        MtrxI   distance{adjacency};
        ssize_t n = distance.rows();
        MtrxI   transshipment(distance.rows(), distance.cols());
        for (ssize_t i = 0; i < n; i++) {
            for (ssize_t j = 0; j < n; j++) {
                transshipment(i, j) = j;
            }
        }
        for (ssize_t i = 0; i < n; i++) {
            distance(i, i)      = 0;
            transshipment(i, i) = 0;
        }

        if (VERBOSE) {
            std::cout << "Distancias iniciales:\n" << distance << "\n\n";
            std::cout << "Transbordos iniciales:\n" << transshipment << "\n\n";
            std::cout << "###################################################################\n";
        }

        for (ssize_t k = 0; k < n; k++) {
            for (ssize_t i = 0; i < n; i++) {
                for (ssize_t j = 0; j < n; j++) {
                    if (distance(i, k) < INF && distance(k, j) < INF &&
                        distance(i, k) + distance(k, j) < distance(i, j)) {
                        distance(i, j)      = distance(i, k) + distance(k, j);
                        transshipment(i, j) = k;
                    }
                }
            }
            if (VERBOSE) {
                std::cout << "k: " << k << '\n';
                std::cout << "Distancias:\n" << distance << "\n\n";
                std::cout << "Transbordos:\n" << transshipment << "\n\n";
                std::cout
                    << "###################################################################\n";
            }
        }
        // RESULTADOS
        std::cout << "Distancias:\n" << distance << "\n\n";
        std::cout << "Transbordos:\n" << transshipment << "\n\n";
        for (ssize_t i = 0; i < adjacency.rows(); i++) {
            std::cout << '[' << i << "] = " << get_node_name(i) << '\n';
        }
        ////////////
    }

    ssize_t GraphProblem::bfs(MtrxI const& capacity) {
        fill(parent.begin(), parent.end(), -1);
        parent[source] = -2;
        std::queue<std::pair<ssize_t, ssize_t>> queue;
        queue.push({source, INF});

        while (!queue.empty()) {
            ssize_t cur  = queue.front().first;
            ssize_t flow = queue.front().second;
            queue.pop();

            for (ssize_t next : adj_list[cur]) {
                if (parent[next] == -1 && capacity(cur, next)) {
                    parent[next] = cur;
                    int new_flow = std::min(flow, capacity(cur, next));
                    if (next == destination) {
                        return new_flow;
                    }
                    queue.push({next, new_flow});
                }
            }
        }

        return 0;
    }

    void GraphProblem::max_flow() {
        parent.resize(adjacency.rows());
        adj_list.resize(adjacency.rows());
        fill(adj_list.begin(), adj_list.end(), std::vector<ssize_t>());

        MtrxI capacity(adjacency);
        for (ssize_t i = 0; i < capacity.rows(); i++) {
            for (ssize_t j = 0; j < capacity.cols(); j++) {
                if (capacity(i, j) != INF && i != j) {
                    adj_list[i].push_back(j);
                    adj_list[j].push_back(i);
                } else if (capacity(i, j) == INF) {
                    capacity(i, j) = 0;
                }
            }
        }

        ssize_t flow = 0;
        for (ssize_t new_flow = bfs(capacity); new_flow; new_flow = bfs(capacity)) {
            std::vector<ssize_t> path;

            flow += new_flow;
            ssize_t current = destination;
            path.emplace(path.begin(), destination);
            while (current != source) {
                ssize_t previous = parent[current];
                path.emplace(path.begin(), previous);
                capacity(previous, current) -= new_flow;
                capacity(current, previous) += new_flow;
                current = previous;
            }
            if (VERBOSE) {
                std::cout << "Flow " << new_flow << ": ";
                for (size_t i = 0; i < path.size(); i++) {
                    if (i + 1 == path.size()) {
                        std::cout << get_node_name(path[i]) << '\n';
                    } else {
                        std::cout << get_node_name(path[i]) << "->";
                    }
                }
            }
        }
        // RESULTADOS
        std::cout << "Max flow: " << flow << '\n';
        /////////////
    }

    void GraphProblem::solve() {
        switch (type) {
            case TRANSPORTATION: {
                transportation();
            } break;

            case ASSIGNMENT: {
                hungarian_method();
            } break;

            case MIN_SPANNING_TREE: {
                min_spanning_tree();
            } break;

            case SHORTEST_PATH_DIJKSTRA: {
                dijkstra_algorithm();
            } break;

            case SHORTEST_PATH_FLOYD: {
                floyd_warshall_algorithm();
            } break;

            case MAX_FLOW: {
                max_flow();
            } break;

            default:
                break;
        }
    }

}  // namespace optimization