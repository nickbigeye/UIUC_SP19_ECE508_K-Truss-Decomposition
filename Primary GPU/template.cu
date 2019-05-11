#include <cstdio>
#include <cstdlib>
#include <stdio.h>

#include "template.hu"

// neighbor-set-intersection-based triangle counting kernel

__device__ int triangle_counts_set(
        const int *edgeDst,
        const int *rowPtr,
        int *W,
        const int u, const int v) {
    int u_start = rowPtr[u], u_end = rowPtr[u + 1];
    int v_start = rowPtr[v], v_end = rowPtr[v + 1];
    int w1 = edgeDst[u_start];
    int w2 = edgeDst[v_start];
    int tc = 0;
    while (u_start < u_end && v_start < v_end) {
        if (w1 == -1 || w1 < w2) {
            w1 = edgeDst[++u_start];
        }
        if (w2 == -1 || w1 > w2) {
            w2 = edgeDst[++v_start];
        }
        if (w1 != -1 && w2 != -1 && w1 == w2) {
            W[tc] = w1;
            tc += 1;
            w1 = edgeDst[++u_start];
            w2 = edgeDst[++v_start];
        }
    }
    return tc;
}

__global__ void markAll(int *affected, int currEdges, int flag) {
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < currEdges) {
        affected[index] = flag;
    }
}

__global__ void checkAffectedEdges(
        int * affected_data,
        int affected_data_length,
        int * edgeSrc,
        int * edgeDst,
        int * to_delete_data,
        int * e_aff_data,
        int * rowPtr,
        int k
        ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < affected_data_length){
        if (e_aff_data[i] == 1 && edgeSrc[i] != -1 && edgeDst[i] != -1) {
            int src = edgeSrc[i];
            int dst = edgeDst[i];
            int tcSet[200];
            int tc = triangle_counts_set(edgeDst, rowPtr, tcSet, src, dst);
            if (tc < k - 2) {
                to_delete_data[i] = 1;

                for (int l = 0; l < tc; ++l) {
                    int w = tcSet[l];
                    // check if (u, w) is in edge
                    for (int j = rowPtr[src]; j < rowPtr[src + 1]; ++j) {
                        if (edgeDst[j] == w) {
                            affected_data[j] = 1;
                        }
                    }
                    // check if (v, w) is in edge
                    for (int j = rowPtr[dst]; j < rowPtr[dst + 1]; ++j) {
                        if (edgeDst[j] == w) {
                            affected_data[j] = 1;
                        }
                    }
                    // check if (w, u) or (w, v) is in edge
                    for (int j = rowPtr[w]; j < rowPtr[w + 1]; ++j) {
                        if (edgeDst[j] == src || edgeDst[j] == dst) {
                            affected_data[j] = 1;
                        }
                    }
                }
            }
        }
    }
}

int k_truss(
        int *edgeSrc,         //!< node ids for edge srcs
        int *edgeDst,         //!< node ids for edge dsts
        int *rowPtr,          //!< source node offsets in edgeDst
        int numEdges,
        int numRow
) {

    int k = 3;

    int * edgeSrcDevice;
    int * edgeDstDevice;
    int * rowPtrDevice;

    cudaMalloc((void**) &edgeSrcDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &edgeDstDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &rowPtrDevice, numRow * sizeof(int));

    pangolin::Vector<int> affected; // a numEdges long vector containing whether the ith edge is affected
    pangolin::Vector<int> to_delete;
    pangolin::Vector<int> e_aff;

    // initialize all data
    for (int i = 0; i < numEdges; i++) {
        affected.push_back(0);
        to_delete.push_back(0);
        e_aff.push_back(0);
    }

    while (true) {

        // mark all edge to be affected and un-deleted if valid edges

        for (int i = 0; i < numEdges; i++) {
            if (edgeSrc[i] != -1){
                affected.data()[i] = 1;
                to_delete.data()[i] = 0;
                e_aff.data()[i] = 0;
            }
        }

        while (true) {

            bool flag = true;

            cudaMemcpy(edgeSrcDevice, edgeSrc, numEdges * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(edgeDstDevice, edgeDst, numEdges * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(rowPtrDevice, rowPtr, numRow * sizeof(int), cudaMemcpyHostToDevice);

            for (int i = 0; i < numEdges; ++i) {
                e_aff.data()[i] = 0; // first set the initial value of e_aff to be zero
                if (affected.data()[i] == 1) {
                    e_aff.data()[i] = 1; // then update it to be 1 if it's affected
                    flag = false;        // and set break flag to be false
                }
            }

            if (flag) {
                break;
            }

            // mark all e as "not affected"
            dim3 dimBlock(512, 1, 1);
            dim3 dimGrid(ceil(affected.size()/512.0), 1, 1);
            markAll<<<dimGrid, dimBlock>>>(affected.data(), affected.size(), -1);
            cudaDeviceSynchronize();

            // check all affected edges
            checkAffectedEdges<<<dimGrid, dimBlock>>>(
                    affected.data(), affected.size(), edgeSrcDevice, edgeDstDevice,
                    to_delete.data(), e_aff.data(),
                    rowPtrDevice, k);
            cudaDeviceSynchronize();

            // short update
            for (int i = 0; i < numEdges; ++i) {
                if (to_delete.data()[i] == 1){
                    edgeSrc[i] = -1;
                    edgeDst[i] = -1;
                }
            }
        }

        int currEdges = 0;
        for (int i = 0; i < numEdges; ++i){
            if (edgeSrc[i] != -1){
                currEdges += 1;
            }
        }

        if (currEdges > 0) {
            printf("Order of truss is %d. # of edges remaining is %d.\n", k, currEdges);
            k += 1;
            if (k > 5){
                break;
            }
        } else {
            break;
        }
    }
    cudaFree(edgeSrcDevice);
    cudaFree(edgeDstDevice);
    cudaFree(rowPtrDevice);

    return k - 1;
}

int do_k_truss(pangolin::COOView<int> view) {

    // test case for our own
//    printf("For the first graph:\n");
//    int numEdges = 16;
//    int edgeSrc[16] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};
//    int edgeDst[16] = {1, 2, 0, 2, 3, 4, 0, 1, 3, 4, 1, 2, 4, 1, 2, 3};
//    int rowPtr[6] = {0, 2, 6, 10, 13, 16};
//    int numRow = 6;
//    printf("k-truss order for the first graph is: %d\n", k_truss(edgeSrc, edgeDst, rowPtr, numEdges, numRow));
//
//    printf("For the second graph:\n");
//    int numEdges2 = 24;
//    int edgeSrc2[30] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 7};
//    int edgeDst2[30] = {1, 2, 7, 0, 3, 4, 0, 5, 6, 7, 1, 4, 1, 3, 5, 2, 4, 6, 2, 5, 7, 0, 2, 6};
//    int rowPtr2[10] = {0, 3, 6, 10, 12, 15, 18, 21, 24};
//    int numRow2 = 9;
//    printf("k-truss order for the second graph is %d\n", k_truss(edgeSrc2, edgeDst2, rowPtr2, numEdges2, numRow2));

    int numEdges3 = view.nnz();
    int * edgeSrc3 = const_cast<int*>(view.row_ind());
    int * edgeDst3 = const_cast<int*>(view.col_ind());
    int * rowPtr3 = const_cast<int*>(view.row_ptr());
    int numRow3 = view.num_rows();

    printf("k-truss order for the input graph is %d\n", k_truss(edgeSrc3, edgeDst3, rowPtr3, numEdges3, numRow3));
    return 0;
}
