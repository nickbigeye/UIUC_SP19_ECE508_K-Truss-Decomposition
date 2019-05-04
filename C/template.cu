#include <cstdio>
#include <cstdlib>
#include <stdio.h>

#include "template.hu"

// neighbor-set-intersection-based triangle counting kernel

__device__ int triangle_counts_set(
        const int *edgeSrc,
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
        int k,
        int * triangleCounts
        ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < affected_data_length){
        if (e_aff_data[i] == 1) {
            int src = edgeSrc[i];
            int dst = edgeDst[i];
            int tcSet[100];
            int tc = triangle_counts_set(edgeSrc, edgeDst, rowPtr, tcSet, src, dst);
            triangleCounts[i] = tc;
            if (tc < k - 2) {
                to_delete_data[i] = 1;
                for (int j = 0; j < affected_data_length; ++j) {
                    if (edgeSrc[j] == src) {
                        for (int l = 0; l < tc; ++l) {
                            if (edgeDst[j] == tcSet[l]) {
                                affected_data[j] = 1;
                            }
                        }
                    } else if (edgeDst[j] == src) {
                        for (int l = 0; l < tc; ++l) {
                            if (edgeSrc[j] == tcSet[l]) {
                                affected_data[j] = 1;
                            }
                        }
                    } else if (edgeSrc[j] == dst) {
                        for (int l = 0; l < tc; ++l) {
                            if (edgeDst[j] == tcSet[l]) {
                                affected_data[j] = 1;
                            }
                        }
                    } else if (edgeDst[j] == dst) {
                        for (int l = 0; l < tc; ++l) {
                            if (edgeSrc[j] == tcSet[l]) {
                                affected_data[j] = 1;
                            }
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

    int currEdges = numEdges;
    int currRow = numRow;

    while (true) {
        pangolin::Vector<int> affected; // a numEdges long vector containing whether the ith edge is affected
        pangolin::Vector<int> to_delete;
        pangolin::Vector<int> e_aff;

        // mark all e in edges as "affected" and un-deleted
        for (int i = 0; i < currEdges; i++) {
            affected.push_back(1);
            to_delete.push_back(0);
            e_aff.push_back(0);
        }

        while (true) {
            pangolin::Vector<int> triangleCounts;
            int * edgeSrcDevice;
            int * edgeDstDevice;
            int * rowPtrDevice;

            cudaMalloc((void**) &edgeSrcDevice, currEdges * sizeof(int));
            cudaMalloc((void**) &edgeDstDevice, currEdges * sizeof(int));
            cudaMalloc((void**) &rowPtrDevice, currRow * sizeof(int));

            cudaMemcpy(edgeSrcDevice, edgeSrc, currEdges * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(edgeDstDevice, edgeDst, currEdges * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(rowPtrDevice, rowPtr, currRow * sizeof(int), cudaMemcpyHostToDevice);

            bool flag = true;

            for (int i = 0; i < currEdges; ++i) {
                e_aff.data()[i] = 0; // first set the initial value of e_aff to be zero
                if (affected.data()[i] == 1) {
                    e_aff.data()[i] = 1; // then update it to be 1 if it's affected
                    flag = false;        // and set break flag to be false
                }
                triangleCounts.push_back(0);
            }

            if (flag) {
                break;
            }

            // mark all e as "not affected"
            dim3 dimBlock(256, 1, 1);
            dim3 dimGrid(ceil(affected.size()/256.0), 1, 1);
            markAll<<<dimGrid, dimBlock>>>(affected.data(), affected.size(), -1);
            cudaDeviceSynchronize();
//            printf("all elements marked unaffected\n");

//            printf("check for src nodes:\n");
//            for (int i = 0; i < currEdges; ++i){
//                printf("%d, ", edgeSrc[i]);
//            }
//            printf("\n");
//
//            printf("check for dst nodes:\n");
//            for (int i = 0; i < currEdges; ++i){
//                printf("%d, ", edgeDst[i]);
//            }
//            printf("\n");
//
//            printf("check for e_aff data:\n");
//            for (int i=0; i < currEdges; ++i){
//                printf("%d, ", e_aff[i]);
//            }
//            printf("\n");

            // check all affected edges
            checkAffectedEdges<<<dimGrid, dimBlock>>>(
                    affected.data(), affected.size(), edgeSrcDevice, edgeDstDevice,
                    to_delete.data(), e_aff.data(),
                    rowPtrDevice, k, triangleCounts.data());
            cudaDeviceSynchronize();

//            // check for to delete data
//            printf("check for edges to be deleted:\n");
//            for (int i = 0; i < currEdges; ++i){
//                if (to_delete.data()[i] == 1){
//                    printf("%d, ", i);
//                }
//            }
//            printf("\n");
//
//            // check for affected data
//            printf("check for edges to be affected:\n");
//            for (int i = 0; i < currEdges; ++i){
//                if (affected.data()[i] == 1){
//                    printf("%d, ", i);
//                }
//            }
//            printf("\n");
//
//            // check for triangle data
//            printf("check for triangle data:\n");
//            for (int i = 0; i < currEdges; ++i){
//                printf("%d, ", triangleCounts.data()[i]);
//            }
//            printf("\n");

            // short update
            for (int i = 0; i < currEdges; ++i) {
                if (to_delete.data()[i] == 1){
                    edgeSrc[i] = -1;
                    edgeDst[i] = -1;
                }
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
        currRow = newPtr.size();

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
            if (k > 5){
                break;
            }
        } else {
            break;
        }
    }

    return k - 1;
}

int do_k_truss() {

    // test case for our own
    printf("For the first graph:\n");
    int numEdges = 16;
    int edgeSrc[16] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};
    int edgeDst[16] = {1, 2, 0, 2, 3, 4, 0, 1, 3, 4, 1, 2, 4, 1, 2, 3};
    int rowPtr[6] = {0, 2, 6, 10, 13, 16};
    int numRow = 6;
    printf("k-truss order for the first graph is: %d\n", k_truss(edgeSrc, edgeDst, rowPtr, numEdges, numRow));

    printf("For the second graph:\n");
    int numEdges2 = 30;
    int edgeSrc2[30] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8};
    int edgeDst2[30] = {1, 2, 7, 0, 3, 4, 0, 5, 6, 7, 1, 4, 8, 1, 3, 5, 8, 2, 4, 6, 2, 5, 7, 8, 0, 2, 6, 3, 4, 6};
    int rowPtr2[10] = {0, 3, 6, 10, 13, 17, 20, 24, 27, 30};
    int numRow2 = 10;
    printf("k-truss order for the second graph is %d\n", k_truss(edgeSrc2, edgeDst2, rowPtr2, numEdges2, numRow2));
//
//    int numEdges3 = view.nnz();
//    int * edgeSrc3 = const_cast<int*>(view.row_ind());
//    int * edgeDst3 = const_cast<int*>(view.col_ind());
//    int * rowPtr3 = const_cast<int*>(view.row_ptr());
//
//    printf("k-truss order for the input graph is %d\n", k_truss(edgeSrc3, edgeDst3, rowPtr3, numEdges3));
    return 0;
}
