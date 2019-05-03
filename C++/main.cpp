#include <iostream>
#include <list>
#include <iterator>
#include <utility>
#include <unordered_map>
#include <map>
#include <tuple>

using namespace std;

list<pair<unsigned int, unsigned int>> get_edges(list<list<unsigned int>> &l_graph);

unordered_map<unsigned int, list<unsigned int>> get_neighbour(list<list<unsigned int>> &l_graph);

tuple<unsigned int, list<unsigned int>> triangle_count(unordered_map<unsigned int, list<unsigned int>> &neighbour, unsigned int &u, unsigned int &v);

void k_truss(list<list<unsigned int>> &l_graph);

template <class T>
bool in(list<T> l, T x) {
    return (find(l.begin(), l.end(), x) != l.end());
}

list<pair<unsigned int, unsigned int>> get_edges(list<list<unsigned int>> &l_graph){
    list<pair<unsigned int, unsigned int>> edges;

    unsigned int i = 0;
    for (auto Row : l_graph) {
        unsigned int j = 0;
        for (auto col : Row) {
            unsigned int connected = col;
            if (connected && i != j) edges.emplace_back(make_pair(i, j));
            j++;
        }
        i++;
    }
    return edges;
}

unordered_map<unsigned int, list<unsigned int>> get_neighbour(list<list<unsigned int>> &l_graph) {
    unordered_map<unsigned int, list<unsigned int>> neighbour;

    unsigned int i = 0;
    for (auto Row : l_graph) {
        unsigned int j = 0;
        for (auto col : Row) {
            unsigned int connected = col;
            if (connected && i != j) {
                neighbour[i].push_back(j);
            }
            j++;
        }
        i++;
    }
    return neighbour;
}

tuple<unsigned int, list<unsigned int>> triangle_count(unordered_map<unsigned int, list<unsigned int>> &neighbour, unsigned int &u, unsigned int &v) {
    unsigned int ptr1 = 0, ptr2 = 0, res = 0;
    list<unsigned int> res_set;

    while (ptr1 < neighbour[u].size() && ptr2 < neighbour[v].size()) {
        unsigned int u_neighbour = *next(neighbour[u].begin(), ptr1);
        unsigned int v_neighbour = *next(neighbour[v].begin(), ptr2);

        if (u_neighbour == v_neighbour) {
            if (!in(res_set, u_neighbour)) res_set.push_back(u_neighbour);
            res++;
            ptr1++;
            ptr2++;
        }
        else if (u_neighbour < v_neighbour) ptr1++;
        else ptr2++;
    }
    return make_tuple(res, res_set);
}

void k_truss(list<list<unsigned int>> &l_graph) {
    int k = 3;
    map<pair<unsigned int, unsigned int>, bool> affected;
    map<pair<unsigned int, unsigned int>, bool> deleted;


    /*
    list<pair<unsigned int, unsigned int>> edges = get_edges(l_graph);
    cout << "----------print edges\n";
    for(auto e : edges) {
        cout << e.first << ' ' << e.second << endl;
    }
    cout << "---------print affected\n";
    for (auto e : edges) {
        affected[e] = true;
    }
    for (auto E : affected) {
        cout << '(' << E.first.first << ", " << E.first.second << "): " << E.second << endl;
    }

    cout << "----------print neighbour\n";
    unordered_map<unsigned int, list<unsigned int>> neighbour = get_neighbour(l_graph);
    for (auto s : neighbour) {
        cout << s.first << ": [";
        for (auto n : s.second) {
            cout << n << ", ";
        }
        cout << ']' << endl;
    }

    cout << "----------print graph\n";
    for (auto Row : l_graph) {
        for (auto col : Row) {
            cout << col << ' ';
        }
        cout << endl;
    }

    cout << "----------print e_affected\n";
    list<pair<unsigned int, unsigned int>> e_affected;
    for (auto e : edges) {
        if (affected.count(e) != 0 && e.first < e.second) e_affected.push_back(e);
    }
    for (auto e_aff : e_affected) {
        cout << '(' << e_aff.first << ", " << e_aff.second << ')' << endl;
    }

    cout << "----------print tc_tuple\n";
    for (auto e_aff : e_affected) {
        tuple<unsigned int, list<unsigned int>> tc_tuple = triangle_count(neighbour, e_aff.first, e_aff.second);
        unsigned int tc = get<0>(tc_tuple);
        list<unsigned int> tc_set = get<1>(tc_tuple);
        cout << "for (" << e_aff.first << ", " << e_aff.second << "):\n";
        cout << "tc: " << tc << endl << "tc_set: ";
        for (auto intersect : tc_set) {
            cout << intersect << ", ";
        }
        cout << endl << "---\n";
    }*/


    while (true) {
        list<pair<unsigned int, unsigned int>> edges = get_edges(l_graph);
        unordered_map<unsigned int, list<unsigned int>> neighbour = get_neighbour(l_graph);

        cout << "----------print edges\n";
        for(auto e : edges) {
            cout << e.first << ' ' << e.second << endl;
        }

        // mark all e in edges as "affected"
        for (auto e : edges) {
            affected[e] = true;
        }

        while (true) {
            list<pair<unsigned int, unsigned int>> e_affected;
            for (auto e : edges) {
                if (affected.count(e) != 0 && e.first < e.second) e_affected.push_back(e);
            }

            if (e_affected.empty()) break;

            // mark all e in edges as "unaffected"
            affected.clear();

            for (auto e_aff : e_affected) {
                unsigned int u = e_aff.first;
                unsigned int v = e_aff.second;

                tuple<unsigned int, list<unsigned int>> tc_tuple = triangle_count(neighbour, u, v);
                unsigned int tc = get<0>(tc_tuple);
                list<unsigned int> tc_set = get<1>(tc_tuple);

                if (tc < k - 2) {
                    // stage 1, mark e_forward and e_reverse as "delete"
                    pair<unsigned int, unsigned int> e_forward = make_pair(u, v);
                    pair<unsigned int, unsigned int> e_reverse = make_pair(v, u);
                    deleted[e_forward] = true;
                    deleted[e_reverse] = true;

                    // stage 2, mark 4 edges as "affected"
                    list<unsigned int> w = tc_set;
                    if (!w.empty()) {
                        for (auto w_item : w) {
                            affected[make_pair(u, w_item)] = true;
                            affected[make_pair(v, w_item)] = true;
                            if (in(e_affected, make_pair(w_item, u))) affected[make_pair(w_item, u)] = true;
                            if (in(e_affected, make_pair(w_item, v))) affected[make_pair(w_item, v)] = true;
                        }
                    }
                }
            }

            list<pair<unsigned int, unsigned int>> new_edges;

            for (auto e : edges) {
                unsigned int u = e.first;
                unsigned int v = e.second;

                if (deleted.count(e) == 0) new_edges.push_back(e);
                else {
                   *next(next(l_graph.begin(), u)->begin(), v) = 0;
                   *next(next(l_graph.begin(), v)->begin(), u) = 0;
                }
            }

            edges = new_edges; // To be tested
            neighbour = get_neighbour(l_graph);
        }

        if (!edges.empty()) {
            cout << "*************Result************\n";
            cout << "-----k: " << k << endl;
            for(auto e : edges) {
                cout << e.first << "--" << e.second << endl;
            }
            k++;
        }
        else break;
    }
}




int main() {
    int N_node = 5;

    unsigned int graph[5][5] = {
            {1, 1, 1, 0, 0},
            {1, 1, 1, 1, 1},
            {1, 1, 1, 1, 1},
            {0, 1, 1, 1, 1},
            {0, 1, 1, 1, 1}
    };

    list<list<unsigned int>> l_graph;

    for(int i = 0; i < N_node; i++) {
        list<unsigned int> row;
        for (int j = 0; j < N_node; j++) {
            row.push_back(graph[i][j]);
        }
        l_graph.push_back(row);
    }

    k_truss(l_graph);

    return 0;
}