#include <cstdio>
#include <cstdlib>
#include <stdio.h>

#include "template.hu"

// neighbor-set-intersection-based triangle counting kernel

pangolin::Vector<int> triangle_counts(
        int *edgeDst,
        int *rowPtr,
        int u, int v) {
    int u_start = rowPtr[u], u_end = rowPtr[u + 1];
    int v_start = rowPtr[v], v_end = rowPtr[v + 1];
    int w1 = edgeDst[u_start];
    int w2 = edgeDst[v_start];
    pangolin::Vector<int> W;
    while (u_start < u_end and v_start < v_end) {
        if (w1 == -1 || w1 < w2) {
            w1 = edgeDst[++u_start];
        }
        if (w2 == -1 || w1 > w2) {
            w2 = edgeDst[++v_start];
        }
        if (w1 != -1 && w2 != -1 && w1 == w2) {
            W.push_back(w1);
            w1 = edgeDst[++u_start];
            w2 = edgeDst[++v_start];
        }
    }
    return W;
}

int k_truss(
        int *edgeSrc,         //!< node ids for edge srcs
        int *edgeDst,         //!< node ids for edge dsts
        int *rowPtr,          //!< source node offsets in edgeDst
        int numEdges                  //!< how many edges to count triangles for
) {

    int k = 3;

    // initialize the affected vector
    pangolin::Vector<int> affected;
    for (int i = 0; i < numEdges; ++i){
        affected.push_back(0);
    }

    while (true) {

        // mark all e in edges as "affected"
        for (int i = 0; i < numEdges; i++) {
            if (edgeSrc[i] != -1){
                affected.data()[i] = 1;
            }
        }

        while (true) {
            pangolin::Vector<int> e_aff; // contains the index of affected edges
            pangolin::Vector<int> to_delete; // contains the index of to-be-deleted edges

            for (int i = 0; i < numEdges; ++i) {
                if (affected.data()[i] == 1) {
                    e_aff.push_back(i);
                }
            }

            if (e_aff.size() == 0) {
                break;
            }

            //mark all e as "not affected"
            for (int i = 0; i < numEdges; ++i) {
                affected.data()[i] = -1;
            }

            for (int i = 0; i < e_aff.size(); ++i) {
                int src = edgeSrc[e_aff.data()[i]];
                int dst = edgeDst[e_aff.data()[i]];
                if (src != -1 && dst != -1){
                    pangolin::Vector<int> tcSet = triangle_counts(edgeDst, rowPtr, src, dst);
                    int tc = tcSet.size();
                    if (tc < k - 2) {
                        // stage 1
                        to_delete.push_back(e_aff.data()[i]);

                        // stage 2
                        for (int l = 0; l < tc; ++l){
                            int w = tcSet.data()[l];
                            // check if (u, w) is in edge
                            for (int j = rowPtr[src]; j < rowPtr[src+1]; ++j){
                                if (edgeDst[j] == w){
                                    affected.data()[j] = 1;
                                }
                            }
                            // check if (v, w) is in edge
                            for (int j = rowPtr[dst]; j < rowPtr[dst+1]; ++j){
                                if (edgeDst[j] == w){
                                    affected.data()[j] = 1;
                                }
                            }
                            // check if (w, u) or (w, v) is in edge
                            for (int j = rowPtr[w]; j < rowPtr[w+1]; ++j) {
                                if (edgeDst[j] == src || edgeDst[j] == dst){
                                    affected.data()[j] = 1;
                                }
                            }
                        }
                    }
                }
            }
            // short update
            for (int i = 0; i < to_delete.size(); ++i) {
                edgeSrc[to_delete.data()[i]] = -1;
                edgeDst[to_delete.data()[i]] = -1;
            }
        }

        // no long update any more
        int currEdges = 0;
        for (int i = 0; i < numEdges; ++i){
            if (edgeSrc[i] != -1){
                currEdges += 1;
            }
        }

        if (currEdges > 0) {
            printf("Order of truss is %d. Remaining # of edges is %d\n", k, currEdges);
            k += 1;
        } else {
            break;
        }
    }

    return k - 1;
}

int do_k_truss(pangolin::COOView<int> view) {
// REQUIRED

//@@ create a pangolin::Vector (uint64_t) to hold per-edge triangle counts
// Pangolin is backed by CUDA so you do not need to explicitly copy data between host and device.
// You may find pangolin::Vector::data() function useful to get a pointer for your kernel to use.
//        uint64_t numEdges = view.nnz();
//        const uint32_t * edgeSrc = view.row_ind();
//        const uint32_t * edgeDst = view.col_ind();
//        const uint32_t * rowPtr  = view.row_ptr();

// test case for our own
    printf("For the first graph:\n");
    int numEdges = 16;
    int edgeSrc[16] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};
    int edgeDst[16] = {1, 2, 0, 2, 3, 4, 0, 1, 3, 4, 1, 2, 4, 1, 2, 3};
    int rowPtr[6] = {0, 2, 6, 10, 13, 16};
    printf("result for the first graph is: %d\n", k_truss(edgeSrc, edgeDst, rowPtr, numEdges));

    printf("For the second graph:\n");
    int numEdges2 = 30;
    int edgeSrc2[30] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8};
    int edgeDst2[30] = {1, 2, 7, 0, 3, 4, 0, 5, 6, 7, 1, 4, 8, 1, 3, 5, 8, 2, 4, 6, 2, 5, 7, 8, 0, 2, 6, 3, 4, 6};
    int rowPtr2[10] = {0, 3, 6, 10, 13, 17, 20, 24, 27, 30};
    printf("result for the second graph is %d\n", k_truss(edgeSrc2, edgeDst2, rowPtr2, numEdges2));

    printf("For the real input graph - California Road Network:\n");
    int numEdges3 = view.nnz();
    int * edgeSrc3 = const_cast<int*>(view.row_ind());
    int * edgeDst3 = const_cast<int*>(view.col_ind());
    int * rowPtr3 = const_cast<int*>(view.row_ptr());
    int numRow3 = view.num_rows();
    printf("number of edges: %d\nnumber of nodes: %d\n", numEdges3, numRow3);
    printf("k-truss order for the input graph is %d\n", k_truss(edgeSrc3, edgeDst3, rowPtr3, numEdges3));

    return 0;

}