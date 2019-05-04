#include <cstdio>
#include <cstdlib>
#include <stdio.h>
#include <vector>

#include "template.hu"

// neighbor-set-intersection-based triangle counting kernel

__device__ int triangle_counts(
        int *edgeSrc,
        int *edgeDst,
        int *rowPtr,
        int u, int v) {
    int u_start = rowPtr[u], u_end = rowPtr[u + 1];
    int v_start = rowPtr[v], v_end = rowPtr[v + 1];
    int w1 = edgeDst[u_start];
    int w2 = edgeDst[v_start];
    int n_tcSet = 0;
    while (u_start < u_end and v_start < v_end) {
        if (w1 == -1 || w1 < w2) {
            w1 = edgeDst[++u_start];
        }
        if (w2 == -1 || w1 > w2) {
            w2 = edgeDst[++v_start];
        }
        if (w1 != -1 && w2 != -1 && w1 == w2) {
//            tcSet.push_back(w1);
            w1 = edgeDst[++u_start];
            w2 = edgeDst[++v_start];
            n_tcSet++;
        }
    }
    return n_tcSet;
}

__device__ void get_tcSet(
        int *edgeSrc,
        int *edgeDst,
        int *rowPtr,
        int u, int v,
        int *tcSet) {
    int u_start = rowPtr[u], u_end = rowPtr[u + 1];
    int v_start = rowPtr[v], v_end = rowPtr[v + 1];
    int w1 = edgeDst[u_start];
    int w2 = edgeDst[v_start];
    int n_tcSet = 0;
    while (u_start < u_end and v_start < v_end) {
        if (w1 == -1 || w1 < w2) {
            w1 = edgeDst[++u_start];
        }
        if (w2 == -1 || w1 > w2) {
            w2 = edgeDst[++v_start];
        }
        if (w1 != -1 && w2 != -1 && w1 == w2) {
//            tcSet.push_back(w1);
            tcSet[n_tcSet] = w1;
            w1 = edgeDst[++u_start];
            w2 = edgeDst[++v_start];
            n_tcSet++;
        }
    }
}

__global__ static void kernel_affected_mark(pangolin::Vector<int> &e_aff,
                                            pangolin::Vector<int> &affected,
                                            pangolin::Vector<int> &to_delete,
                                            int *edgeSrc,
                                            int *edgeDst,
                                            int *rowPtr,
                                            int k,
                                            int currEdges) {
    unsigned int edge_index = blockIdx.x * blockDim.x + threadIdx.x;
    if (edge_index < e_aff.size()) {
        int src = edgeSrc[e_aff.data()[edge_index]];
        int dst = edgeDst[e_aff.data()[edge_index]];
        int n_tcSet = triangle_counts(edgeSrc, edgeDst, rowPtr, src, dst);
        printf("------k: %d, n_tcSet: %d\n", k, n_tcSet);
        int * tcSet = (int *)malloc(n_tcSet * sizeof(int));
        get_tcSet(edgeSrc, edgeDst, rowPtr, src, dst, tcSet);
        for (int i = 0; i < n_tcSet; i++) {
            printf("%d: %d\n", i, tcSet[i]);
        }
        int tc = n_tcSet;
        if (tc < k - 2) {
            // stage 1
            to_delete.push_back(e_aff.data()[edge_index]);

            // stage 2
            for (int j = 0; j < currEdges; ++j) {
                if (edgeSrc[j] == src) {
                    for (int k = 0; k < n_tcSet; ++k) {
                        if (edgeDst[j] == tcSet[k]) {
                            affected[j] = 1;
                        }
                    }
                } else if (edgeDst[j] == src) {
                    for (int k = 0; k < n_tcSet; ++k) {
                        if (edgeSrc[j] == tcSet[k]) {
                            affected[j] = 1;
                        }
                    }
                } else if (edgeSrc[j] == dst) {
                    for (int k = 0; k < n_tcSet; ++k) {
                        if (edgeDst[j] == tcSet[k]) {
                            affected[j] = 1;
                        }
                    }
                } else if (edgeDst[j] == dst) {
                    for (int k = 0; k < n_tcSet; ++k) {
                        if (edgeSrc[j] == tcSet[k]) {
                            affected[j] = 1;
                        }
                    }
                }
            }
        }
        free(tcSet);
    }
}

int k_truss(
        int *edgeSrc,         //!< node ids for edge srcs
        int *edgeDst,         //!< node ids for edge dsts
        int *rowPtr,          //!< source node offsets in edgeDst
        int numEdges                  //!< how many edges to count triangles for
) {

    int k = 3;

    int currEdges = numEdges;

    while (true) {
        pangolin::Vector<int> affected; // a numEdges long vector containing whether the ith edge is affected

        // mark all e in edges as "affected"
        for (int i = 0; i < currEdges; i++) {
            affected.push_back(1);
        }

        while (true) {
            pangolin::Vector<int> e_aff; // contains the index of affected edges
            pangolin::Vector<int> to_delete; // contains the index of to-be-deleted edges

            for (int i = 0; i < currEdges; ++i) {
                if (affected.data()[i] == 1) {
                    e_aff.push_back(i);
                }
            }

            if (e_aff.size() == 0) {
                break;
            }

            //mark all e as "not affected"
            for (int i = 0; i < currEdges; ++i) {
                affected.data()[i] = -1;
            }

//            for (int i = 0; i < e_aff.size(); ++i) {
//                int src = edgeSrc[e_aff.data()[i]];
//                int dst = edgeDst[e_aff.data()[i]];
//                pangolin::Vector<int> tcSet = triangle_counts(edgeSrc, edgeDst, rowPtr, src, dst);
//                int tc = tcSet.size();
//                if (tc < k - 2) {
//                    // stage 1
//                    to_delete.push_back(e_aff.data()[i]);
//
//                    // stage 2
//                    for (int j = 0; j < currEdges; ++j) {
//                        if (edgeSrc[j] == src) {
//                            for (int k = 0; k < tcSet.size(); ++k) {
//                                if (edgeDst[j] == tcSet.data()[k]) {
//                                    affected[j] = 1;
//                                }
//                            }
//                        } else if (edgeDst[j] == src) {
//                            for (int k = 0; k < tcSet.size(); ++k) {
//                                if (edgeSrc[j] == tcSet.data()[k]) {
//                                    affected[j] = 1;
//                                }
//                            }
//                        } else if (edgeSrc[j] == dst) {
//                            for (int k = 0; k < tcSet.size(); ++k) {
//                                if (edgeDst[j] == tcSet.data()[k]) {
//                                    affected[j] = 1;
//                                }
//                            }
//                        } else if (edgeDst[j] == dst) {
//                            for (int k = 0; k < tcSet.size(); ++k) {
//                                if (edgeSrc[j] == tcSet.data()[k]) {
//                                    affected[j] = 1;
//                                }
//                            }
//                        }
//                    }
//                }
//            }
            dim3 dimBlock(512);
            dim3 dimGrid (ceil(1.0 * e_aff.size() / 512), 1, 1);
            kernel_affected_mark<<<dimGrid, dimBlock>>>(e_aff, affected, to_delete, edgeSrc, edgeDst, rowPtr, k, currEdges);
            cudaDeviceSynchronize();
            // short update
            for (int i = 0; i < to_delete.size(); ++i) {
                edgeSrc[to_delete.data()[i]] = -1;
                edgeDst[to_delete.data()[i]] = -1;
            }
        }

        // long update
        pangolin::Vector<int> newSrc;
        pangolin::Vector<int> newDst;
        pangolin::Vector<int> newPtr;
        for (int i = 0; i < currEdges; ++i) {
            if (edgeSrc[i] != -1) {
                newSrc.push_back(edgeSrc[i]);
            }
            if (edgeDst[i] != -1) {
                newDst.push_back(edgeDst[i]);
            }
        }

        newPtr.push_back(0);
        for (int i = 1; i < newSrc.size(); ++i) {
            if (newSrc.data()[i] != newSrc.data()[i - 1]) {
                newPtr.push_back(i);
            }
        }
        newPtr.push_back(newSrc.size());

        // have deep copy every thing to edgeSrc, edgeDst and rowPtr
        edgeSrc = (int *) malloc(newSrc.size() * sizeof(int));
        for (int i = 0; i < newSrc.size(); ++i) {
            edgeSrc[i] = newSrc.data()[i];
        }
        edgeDst = (int *) malloc(newDst.size() * sizeof(int));
        for (int i = 0; i < newDst.size(); ++i) {
            edgeDst[i] = newDst.data()[i];
        }
        rowPtr = (int *) malloc(newPtr.size() * sizeof(int));
        for (int i = 0; i < newPtr.size(); ++i) {
            rowPtr[i] = newPtr.data()[i];
        }
        currEdges = newSrc.size();

        if (currEdges > 0) {
            printf("Order of truss is %d.\n", k);
            printf("src nodes are:\n");
            for (int i = 0; i < newSrc.size(); ++i) {
                printf("%d, ", edgeSrc[i]);
            }
            printf("\n");

            printf("dst nodes are:\n");
            for (int i = 0; i < newDst.size(); ++i) {
                printf("%d, ", edgeDst[i]);
            }
            printf("\n");
            k += 1;
        } else {
            break;
        }
    }

    return k - 1;
}

//__global__ static void kernel_tc(uint64_t *__restrict__ triangleCounts, //!< per-edge triangle counts
//                                 const uint32_t *const edgeSrc,         //!< node ids for edge srcs
//                                 const uint32_t *const edgeDst,         //!< node ids for edge dsts
//                                 const uint32_t *const rowPtr,          //!< source node offsets in edgeDst
//                                 const size_t numEdges                  //!< how many edges to count
//                                                                        //!< triangles for
//) {
//    int index = blockIdx.x * blockDim.x + threadIdx.x;
//
//    if (index < numEdges){
//        uint64_t count = 0;
//
//        // Determine the source and destination node for the edge
//
//        uint32_t src = edgeSrc[index]; // u
//        uint32_t dst = edgeDst[index]; // v
//
//        // Use the row pointer array to determine the start and
//        // end of the neighbor list in the column index array
//
//        uint32_t ptrSrc = rowPtr[src];
//        uint32_t endSrc = rowPtr[src+1];
//
//        uint32_t ptrDst = rowPtr[dst];
//        uint32_t endDst = rowPtr[dst+1];
//
//
//        // Determine how many elements of those two arrays are common
//
//        uint32_t w1 = edgeDst[ptrSrc];
//        uint32_t w2 = edgeDst[ptrDst];
//        while (ptrSrc < endSrc and ptrDst < endDst) {
//            if (w1 < w2) {
//                w1 = edgeDst[++ptrSrc];
//            } else if (w1 > w2) {
//                w2 = edgeDst[++ptrDst];
//            } else if (w1 == w2) {
//                w1 = edgeDst[++ptrSrc];
//                w2 = edgeDst[++ptrDst];
//                count += 1;
//            }
//        }
//        triangleCounts[index] = count;
//    }
//}

uint64_t count_triangles(const pangolin::COOView <uint32_t> view, const int mode) {
    if (mode == 1) {

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

        return 0;
    } else {
        return 0;
    }
}
