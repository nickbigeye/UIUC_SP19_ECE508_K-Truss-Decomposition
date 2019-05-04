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

__global__ void shortUpdate(int *edgeSrc, int *edgeDst, int *to_delete_data, int currEdge){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("Hello from short update!\n");
    if (i < currEdge){
        if (to_delete_data[i] == 1){
            to_delete_data[i] = 0;
            edgeSrc[i] = -1;
            edgeDst[i] = -1;
        }
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
//    printf("From CheckAffectedEdges!\n");
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < affected_data_length){
        if (e_aff_data[i] == 1 && edgeSrc[i] != -1 && edgeDst[i] != -1) {
            int src = edgeSrc[i];
            int dst = edgeDst[i];
            int tcSet[200];
            int tc = triangle_counts_set(edgeSrc, edgeDst, rowPtr, tcSet, src, dst);

            if (tc < k - 2) {
//                printf("%dth edge marked as delete, with tc %d\n", i, tc);
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

__global__ void selectAff(int * e_aff_data, int * affected_data, int currEdges, int * length){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < currEdges){
        e_aff_data[i] = 0; // first set the initial value of e_aff to be zero
        if (affected_data[i] == 1) {
//            printf("%dth edge marked as affected!\n", i);
            e_aff_data[i] = 1; // then update it to be 1 if it's affected
            atomicAdd(length, 1);
        }
    }
}

__global__ void mark(int *edgeSrc, int *affected, int *e_aff, int *to_delete, int numEdges){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("From mark!");
    if (i < numEdges){
        if (edgeSrc[i] != -1){
            affected[i] = 1;
            to_delete[i] = 0;
            e_aff[i] = 0;
        }
    }
}

__global__ void countEdges(int *edgeSrc, int numEdges, int *result){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
//    printf("From countEdges!");
    if (i < numEdges){
        if (edgeSrc[i] != -1){
//            printf("from countEdges: %d, ", i);
            atomicAdd(result, 1);
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

    int * length = (int *) malloc (sizeof(int));
    int * edgeCount = (int *) malloc (sizeof(int));

    int * edgeSrcDevice;
    int * edgeDstDevice;
    int * rowPtrDevice;
    int * lengthDevice;
    int * edgeCountDevice;
    int * affectedDevice;
    int * to_deleteDevice;
    int * e_affDevice;

    cudaMalloc((void**) &edgeSrcDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &edgeDstDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &affectedDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &to_deleteDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &e_affDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &rowPtrDevice, numRow * sizeof(int));
    cudaMalloc((void**) &lengthDevice, 1 * sizeof(int)); // keep track of the length of e_aff;
    cudaMalloc((void**) &edgeCountDevice, 1 * sizeof(int)); // keep track of the length of E after long update;

    cudaMemcpy(edgeSrcDevice, edgeSrc, numEdges * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(edgeDstDevice, edgeDst, numEdges * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(rowPtrDevice,  rowPtr,  numRow * sizeof(int),   cudaMemcpyHostToDevice);


    dim3 dimBlock(512, 1, 1);
    dim3 dimGrid(ceil(numEdges/512.0), 1, 1);

    while (true) {
        // first mark all edge available to be affected and un-deleted
        mark<<<dimGrid, dimBlock>>>(edgeSrcDevice, affectedDevice, e_affDevice, to_deleteDevice, numEdges);
//        printf("Pre-processing all data!\n");

        while (true) {

            // select e_aff out of affected
            length[0] = 0;
            cudaMemcpy(lengthDevice, length, sizeof(int), cudaMemcpyHostToDevice);
            selectAff<<<dimGrid, dimBlock>>>(e_affDevice, affectedDevice, numEdges, lengthDevice);
            cudaDeviceSynchronize();
            cudaMemcpy(length, lengthDevice, sizeof(int), cudaMemcpyDeviceToHost);
//            printf("select aff!\n");

            if (length[0] == 0) {
                break;
            }

            // mark all e as "not affected"
            markAll<<<dimGrid, dimBlock>>>(affectedDevice, numEdges, -1);
            cudaDeviceSynchronize();

//            printf("Mark all as unaffected!\n");

            // check all affected edges
            checkAffectedEdges<<<dimGrid, dimBlock>>>(
                    affectedDevice, numEdges, edgeSrcDevice, edgeDstDevice,
                    to_deleteDevice, e_affDevice, rowPtrDevice, k);
            cudaDeviceSynchronize();

//            printf("Checked all affected edges!\n");

            // short update
            shortUpdate<<<dimGrid, dimBlock>>>(edgeSrcDevice, edgeDstDevice, to_deleteDevice, numEdges);
            cudaDeviceSynchronize();

//            printf("Short update done!\n");
        }

        // long update
        edgeCount[0] = 0;
        cudaMemcpy(edgeCountDevice, edgeCount, 1 * sizeof(int), cudaMemcpyHostToDevice);
        countEdges<<<dimGrid, dimBlock>>>(edgeSrcDevice, numEdges, edgeCountDevice);
        cudaDeviceSynchronize();
//        printf("\n");
        cudaMemcpy(edgeCount, edgeCountDevice, 1 * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(edgeSrc, edgeSrcDevice, numEdges * sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(edgeDst, edgeDstDevice, numEdges * sizeof(int), cudaMemcpyDeviceToHost);

//        printf("Long update done!\n");

//        printf("Src nodes:\n");
//        for (int i = 0; i < numEdges; ++i){
//            printf("%d, ", edgeSrc[i]);
//        }
//
//        printf("\n");
//        printf("Dst nodes:\n");
//        for (int i = 0; i < numEdges; ++i){
//            printf("%d, ", edgeDst[i]);
//        }
//        printf("\n");

        if (edgeCount[0] > 0) {
            printf("Order of truss is %d. Remaining # of edges is %d\n", k, edgeCount[0]);
            k += 1;
            if (k > 10){
                break;
            }
        } else {
            k -= 1;
            break;
        }

    }

    cudaFree(edgeSrcDevice);
    cudaFree(edgeDstDevice);
    cudaFree(rowPtrDevice);
    cudaFree(lengthDevice);
    cudaFree(edgeCountDevice);

    return k;
}

int do_k_truss(pangolin::COOView<int> view) {

//    // test case for our own
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

    printf("For the third graph - California Road Network:\n");
    int numEdges3 = view.nnz();
    int * edgeSrc3 = const_cast<int*>(view.row_ind());
    int * edgeDst3 = const_cast<int*>(view.col_ind());
    int * rowPtr3 = const_cast<int*>(view.row_ptr());
    int numRow3 = view.num_rows();
    printf("number of edges: %d\nnumber of rows: %d\n", numEdges3, numRow3);
    printf("k-truss order for the input graph is %d\n", k_truss(edgeSrc3, edgeDst3, rowPtr3, numEdges3, numRow3));
    return 0;
}
