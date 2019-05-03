#include <cstdio>
#include <cstdlib>
#include <stdio.h>

#include "template.hu"

// neighbor-set-intersection-based triangle counting kernel

pangolin::Vector<int> triangle_counts(
        int *edgeSrc,
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
            w2 = edgeSrc[++v_start];
        }
    }
    return W;
}

int inArr(int value, int *arr, int arrSize) {
    for (int i = 0; i < arrSize; ++i) {
        if (arr[i] == value) {
            return 1;
        }
    }
    return 0;
}

int k_truss(
        int *edgeSrc,         //!< node ids for edge srcs
        int *edgeDst,         //!< node ids for edge dsts
        int *rowPtr,          //!< source node offsets in edgeDst
        int numEdges                  //!< how many edges to count triangles for
) {

    int k = 3;

    while (true) {
        pangolin::Vector<int> affected; // a numEdges long vector containing whether the ith edge is affected

        // mark all e in edges as "affected"
        for (int i = 0; i < numEdges; i++) {
            affected.push_back(1);
        }

        while (true) {
            pangolin::Vector<int> e_aff; // contains the index of affected edges
            pangolin::Vector<int> to_delete; // contains the index of to-be-deleted edges

            for (int i = 0; i < numEdges; ++i) {
                int src = edgeSrc[i];
                int dst = edgeDst[i];
                if (affected.data()[i] == 1 && src < dst) {
                    e_aff.push_back(i);
                }
            }

            if (e_aff.size() == 0) {
                break;
            }

            //mark all e as "not affected"
            for (int i = 0; i < numEdges; ++i){
                affected.data()[i] = -1;
            }

            for (int i = 0; i < e_aff.size(); ++i) {
                int src = edgeSrc[i];
                int dst = edgeDst[i];
                pangolin::Vector<int> tcSet = triangle_counts(edgeSrc, edgeDst, rowPtr, src, dst);
                int tc = tcSet.size();
                if (tc < k - 2) {
                    // stage 1
                    to_delete.push_back(e_aff.data()[i]);

                    // stage 2
                    for (int j = 0; j < tcSet.size(); ++j) {
                        int w_item = tcSet.data()[j];

                        if ((inArr(w_item, e_aff_src.data(), e_aff_src.size()) == 1 &&
                             inArr(src, e_aff_dst.data(), e_aff_dst.size()) == 1) ||
                            (inArr(src, e_aff_src.data(), e_aff_src.size()) == 1 &&
                             inArr(w_item, e_aff_dst.data(), e_aff_dst.size()) == 1)) {
                            affected_src.push_back(w_item);
                            affected_dst.push_back(src);
                            affected_src.push_back(src);
                            affected_dst.push_back(w_item);
                        }

                        if ((inArr(w_item, e_aff_src.data(), e_aff_src.size()) == 1 &&
                             inArr(dst, e_aff_dst.data(), e_aff_dst.size()) == 1) ||
                            (inArr(dst, e_aff_src.data(), e_aff_src.size()) == 1 &&
                             inArr(w_item, e_aff_dst.data(), e_aff_dst.size()) == 1)) {
                            affected_src.push_back(w_item);
                            affected_dst.push_back(dst);
                            affected_src.push_back(dst);
                            affected_dst.push_back(w_item);
                        }
                    }
                }
            }

            // short update
            for (int i = 0; i < to_delete.size(); ++i){
                edgeSrc[to_delete.data()[i]] = -1;
                edgeDst[to_delete.data()[i]] = -1;
            }
        }

        // long update
        pangolin::Vector<int> newSrc;
        pangolin::Vector<int> newDst;
        pangolin::Vector<int> newPtr;
        for (int i = 0; i < numEdges; ++i){
            if (edgeSrc[i] != -1){
                newSrc.push_back(edgeSrc[i]);
            }
            if (edgeDst[i] != -1){
                newDst.push_back(edgeDst[i]);
            }
        }

        newPtr.push_back(0);
        for (int i = 1; i < newSrc.size(); ++i){
            if (newSrc.data()[i] != newSrc.data()[i-1]){
                newPtr.push_back(i);
            }
        }
        newPtr.push_back(newSrc.size());

        edgeSrc = newSrc.data();
        edgeDst = newDst.data();
        rowPtr = newPtr.data();
        numEdges = newSrc.size();

        if (numEdges > 0) {
            k += 1;
        } else {
            break;
        }
    }

    return k;
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
        int numEdges = 8;
        int edgeSrc[16] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};
        int edgeDst[16] = {1, 2, 0, 2, 3, 4, 0, 1, 3, 4, 1, 2, 4, 1, 2, 3};
        int rowPtr[6] = {0, 2, 6, 10, 13, 16};

        return k_truss(edgeSrc, edgeDst, rowPtr, numEdges);
    } else {
        return 0;
    }
}
