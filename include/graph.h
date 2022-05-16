#ifndef GRAPH_H
#define GRAPH_H
#include <Eigen/Dense>
#include <set>

namespace optimization {
    typedef Eigen::Matrix<long long, Eigen::Dynamic, Eigen::Dynamic> MtrxI;
    typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>      MtrxB;
    typedef std::pair<ssize_t, ssize_t>                              pair_idx;
    class GraphProblem {
       public:
        GraphProblem();
        GraphProblem(char const* _name);
        void solve();
        void load_problem(std::ifstream& file);

       private:
        const long long INF = std::numeric_limits<long long>::max() / 2;

        enum GProblem_T {
            NOT_ASSIGNED,
            TRANSPORTATION,
            ASSIGNMENT,
            MIN_SPANNING_TREE,
            SHORTEST_PATH_DIJKSTRA,
            SHORTEST_PATH_FLOYD,
            MAX_FLOW,
            MIN_FLOW,  // no dispobible
        };
        // Cualquier nodo que no sea de suministro o demanda es de tipo default
        enum Node_T { SUPPLY, DEMAND, TRANSHIPMENT, DEFAULT };

        struct Node {
            Node(std::string const& _name = "", ssize_t _index = -1, Node_T _type = DEFAULT)
                : name(_name), index(_index), type(_type) {}
            std::string name;
            ssize_t     index;
            Node_T      type;
        };
        struct Edge {
            ssize_t u;
            ssize_t v;
            ssize_t weight;  // peso de u a v

            bool operator<(Edge const& other) { return weight < other.weight; }
        };
        struct CompareName {
            using is_transparent = std::string;

            bool operator()(Node const& foo_1, Node const& foo_2) const {
                return foo_1.name < foo_2.name;
            }
            bool operator()(std::string const& id, Node const& foo) const { return id < foo.name; }
            bool operator()(Node const& foo, std::string const& id) const { return foo.name < id; }
        };

        void transportation();

        void hungarian_method();

        void min_spanning_tree();

        void dijkstra_algorithm();

        void floyd_warshall_algorithm();  // LISTO

        void max_flow();  // LISTO

        void min_flow();

        void add_dummy_node();

        int is_not_balanced();

        void build_tableau(MtrxI& costs, MtrxI& tableau, std::vector<ssize_t>& supply_relation,
                           std::vector<ssize_t>& demand_relation) const;

        void vogel_method(MtrxI& costs, MtrxI& tableau, MtrxB& basic_variables) const;

        bool bfs_not_found(std::vector<bool> const& supply_available,
                           std::vector<bool> const& demand_available,
                           size_t& total_supply_available, size_t& total_demand_available) const;

        pair_idx is_optimal(MtrxI const& costs, MtrxB const& basic_variables) const;

        void find_loop(pair_idx const& idx, MtrxB const& basic_variables,
                       std::vector<pair_idx>& loop) const;

        void build_hungarian_tableau(MtrxI& costs, std::vector<ssize_t>& supply_relation,
                                     std::vector<ssize_t>& demand_relation) const;

        void subtract_min(MtrxI& costs) const;

        void draw_lines(MtrxI const& costs, MtrxI& lines,
                        std::vector<ssize_t> const& match_col) const;

        ssize_t find_min_lines(MtrxI const& costs, MtrxI& lines) const;

        bool try_kuhn(ssize_t i, std::vector<bool>& used, std::vector<ssize_t>& match_col,
                      MtrxI const& costs) const;

        size_t find_set(ssize_t v);

        void union_sets(ssize_t a, ssize_t b);

        std::vector<ssize_t> restore_path(ssize_t s, ssize_t t,
                                          std::vector<ssize_t> const& p) const;

        std::string get_node_name(ssize_t idx) const;

        ssize_t bfs(MtrxI const& capacity);

        std::string                 name{""};
        GProblem_T                  type{NOT_ASSIGNED};
        std::set<Node, CompareName> nodes;  // Cambiar a unordened_set
        MtrxI                       adjacency;

        // Variables auxiliares para algoritmo de shortest-path
        ssize_t source{-1};
        ssize_t destination{-1};

        // Variable aux para el problema de transporte
        std::vector<long long> b_values;
        long long              total_supply{0};
        long long              total_demand{0};
        size_t                 num_supply_pnts{0};
        size_t                 num_demand_pnts{0};

        // Variables aux para min spanning tree
        std::vector<ssize_t> rank;
        std::vector<ssize_t> parent;

        // Variable aux para max flow
        std::vector<std::vector<ssize_t>> adj_list;
    };
}  // namespace optimization

#endif