import copy


def get_edges(graph):
    edges = list()
    for (i, arr) in enumerate(graph):
        for (j, connected) in enumerate(arr):
            if connected and i != j:
                edges.append((i, j))
    return edges


def get_neighbour(graph):
    neighbour = dict()
    for (i, arr) in enumerate(graph):
        for (j, connected) in enumerate(arr):
            if connected and i != j:
                try:
                    neighbour[i].append(j)
                except KeyError:
                    neighbour[i] = [j]
    return neighbour


def triangle_count(neighbour, u, v):
    ptr1 = 0
    ptr2 = 0
    res = 0
    res_set = list()
    while ptr1 < len(neighbour[u]) and ptr2 < len(neighbour[v]):
        if neighbour[u][ptr1] == neighbour[v][ptr2]:
            if neighbour[u][ptr1] not in res_set:
                res_set.append(neighbour[u][ptr1])
            res += 1
            ptr1 += 1
            ptr2 += 1
        elif neighbour[u][ptr1] < neighbour[v][ptr2]:
            ptr1 += 1
        else:
            ptr2 += 1

    if res != len(res_set):
        print("Not equal!")
    return res, res_set


def k_truss(graph):
    k = 3
    affected = dict()
    deleted = dict()

    while True:

        # long update
        edges = get_edges(graph)
        neighbour = get_neighbour(graph)

        # mark all e in edges as "affected"
        for e in edges:
            affected[e] = True

        while True:
            # the directed edge of affected undirected edges
            e_affected = list()

            for e in edges:
                if e in affected and e[0] < e[1]:
                    e_affected.append(e)

            if len(e_affected) == 0:
                break

            # mark all e in edges as "unaffected"
            affected = dict()

            for e in e_affected:
                tc, tc_set = triangle_count(neighbour, e[0], e[1])
                # print(e, tc)
                if tc < k - 2:
                    # stage 1, mark e_forward and e_reverse as "delete"
                    e_forward = (e[0], e[1])
                    e_reverse = (e[1], e[0])
                    deleted[e_forward] = True
                    deleted[e_reverse] = True

                    # stage 2, mark 4 edges as "affected"
                    w = tc_set
                    for w_item in w:
                        if (w_item, e[0]) in e_affected or (e[0], w_item) in e_affected:
                            affected[(w_item, e[0])] = True
                            affected[(e[0], w_item)] = True
                        if (w_item, e[1]) in e_affected or (e[1], w_item) in e_affected:
                            affected[(w_item, e[1])] = True
                            affected[(e[1], w_item)] = True

            new_edge = list()

            for e in edges:
                if e in deleted:
                    # if labeled deleted, mark as the connected as 0
                    graph[e[0]][e[1]] = 0
                    graph[e[1]][e[0]] = 0
                else:
                    # else, add edge to new edge list
                    new_edge.append(e)

            # short update
            edges = copy.deepcopy(new_edge)
            neighbour = get_neighbour(graph)

        if edges:
            print(k)
            print(edges)
            k += 1
        else:
            break


def main():
    graph = [
        [1, 1, 1, 0, 0],
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],
        [0, 1, 1, 1, 1],
        [0, 1, 1, 1, 1]
    ]

    graph2 = [
        [1, 1, 1, 0, 0, 0, 0, 1, 0],
        [1, 1, 0, 1, 1, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 1, 1, 1, 0],
        [0, 1, 0, 1, 1, 0, 0, 0, 1],
        [0, 1, 0, 1, 1, 1, 0, 0, 1],
        [0, 0, 1, 0, 1, 1, 1, 0, 0],
        [0, 0, 1, 0, 0, 1, 1, 1, 1],
        [1, 0, 1, 0, 0, 0, 1, 1, 0],
        [0, 0, 0, 1, 1, 0, 1, 0, 1]
    ]
    print("For first graph:")
    k_truss(graph)

    print("For second graph:")
    k_truss(graph2)


if __name__ == "__main__":
    main()
