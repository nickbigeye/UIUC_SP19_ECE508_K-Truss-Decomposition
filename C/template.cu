#include <cstdio>
#include <cstdlib>
#include <stdio.h>

#include "template.hu"

// neighbor-set-intersection-based triangle counting kernel

#define BLOCK_SIZE 512

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
//    __global__ void shortUpdate(int *edgeSrc, int *edgeDst, int *to_delete_data, int *edgeFlag, int currEdge){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < currEdge){
        if (to_delete_data[i] == 1){
            to_delete_data[i] = 0;
//            edgeFlag[i] = 0;
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
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < affected_data_length){
        if (e_aff_data[i] == 1 && edgeSrc[i] != -1 && edgeDst[i] != -1) {
            int src = edgeSrc[i]; // u
            int dst = edgeDst[i]; // v
            int tcSet[200];
            int tc = triangle_counts_set(edgeSrc, edgeDst, rowPtr, tcSet, src, dst);
//            printf("Triangle counting for %dth edge is: %d\n.", i, tc);

            if (tc < k - 2) {
                to_delete_data[i] = 1;

                for (int l = 0; l < tc; ++l){
                    int w = tcSet[l];
                    // check if (u, w) is in edge
                    for (int j = rowPtr[src]; j < rowPtr[src+1]; ++j){
                        if (edgeDst[j] == w){
                            affected_data[j] = 1;
                        }
                    }
                    // check if (v, w) is in edge
                    for (int j = rowPtr[dst]; j < rowPtr[dst+1]; ++j){
                        if (edgeDst[j] == w){
                            affected_data[j] = 1;
                        }
                    }
                    // check if (w, u) or (w, v) is in edge
                    for (int j = rowPtr[w]; j < rowPtr[w+1]; ++j) {
                        if (edgeDst[j] == src || edgeDst[j] == dst){
                            affected_data[j] = 1;
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

__global__ void flagInit(int *edgeFlag, int numEdges){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < numEdges) edgeFlag[i] = 1;
}

__global__ void scan(int *X, int *Y, int *S, int numEdges) {
    //@@ INSERT CODE HERE
    __shared__ int XY[2*BLOCK_SIZE]; // block size change
    int i = 2*blockIdx.x*blockDim.x + threadIdx.x;

    if (i < numEdges) XY[threadIdx.x] = X[i];
    if (i+blockDim.x < numEdges) XY[threadIdx.x+blockDim.x] = X[i+blockDim.x];
    for (unsigned int stride = 1; stride <= blockDim.x; stride *= 2) {
        __syncthreads();
        int index = (threadIdx.x+1) * 2* stride -1;
        if (index < 2*blockDim.x) XY[index] += XY[index - stride];
    }
    for (int stride = blockDim.x/2; stride > 0; stride /= 2) {
        __syncthreads();
        int index = (threadIdx.x+1)*stride*2 - 1;
        if(index + stride < 2*blockDim.x) {
            XY[index + stride] += XY[index];
        }
    }
    __syncthreads();
    if (i < numEdges) Y[i+1] = XY[threadIdx.x];
    if (i+blockDim.x < numEdges) Y[i+blockDim.x+1] = XY[threadIdx.x+blockDim.x];

    __syncthreads();
    if (threadIdx.x == blockDim.x - 1) {
        S[blockIdx.x] = XY[2*BLOCK_SIZE - 1];
    }

}

__global__ void scan2(int *X, int *Y, int numEdges) {
    //@@ INSERT CODE HERE
    __shared__ int XY[2*BLOCK_SIZE]; // block size change
    int i = 2*blockIdx.x*blockDim.x + threadIdx.x;

    if (i < numEdges) XY[threadIdx.x] = X[i];
    if (i+blockDim.x < numEdges) XY[threadIdx.x+blockDim.x] = X[i+blockDim.x];
    for (unsigned int stride = 1; stride <= blockDim.x; stride *= 2) {
        __syncthreads();
        int index = (threadIdx.x+1) * 2* stride -1;
        if (index < 2*blockDim.x) XY[index] += XY[index - stride];
    }
    for (int stride = blockDim.x/2; stride > 0; stride /= 2) {
        __syncthreads();
        int index = (threadIdx.x+1)*stride*2 - 1;
        if(index + stride < 2*blockDim.x) {
            XY[index + stride] += XY[index];
        }
    }
    __syncthreads();
    if (i < numEdges) Y[i+1] = XY[threadIdx.x];
    if (i+blockDim.x < numEdges) Y[i+blockDim.x+1] = XY[threadIdx.x+blockDim.x];
}

__global__ void addScanBlockSum(int *S, int *Y, int numEdges) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < numEdges) {
        int temp = Y[i];
        Y[i] = S[blockIdx.x-1] + temp;
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

//    int * length = (int *) malloc (sizeof(int));
//    int * edgeCount = (int *) malloc (sizeof(int));

    pangolin::Vector<int> length;
    length.push_back(0);
    pangolin::Vector<int> edgeCount;
    edgeCount.push_back(0);

    int * edgeSrcDevice;
    int * edgeDstDevice;
    int * rowPtrDevice;
//    int * lengthDevice;
//    int * edgeCountDevice;
    int * affectedDevice;
    int * to_deleteDevice;
    int * e_affDevice;
//    int * edgeFlagDevice;
//    int * edgePreSumDevice;

    cudaMalloc((void**) &edgeSrcDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &edgeDstDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &affectedDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &to_deleteDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &e_affDevice, numEdges * sizeof(int));
    cudaMalloc((void**) &rowPtrDevice, numRow * sizeof(int));
//    cudaMalloc((void**) &edgeFlagDevice, numEdges * sizeof(int));
//    cudaMalloc((void**) &edgePreSumDevice, numEdges * sizeof(int));

//    cudaMalloc((void**) &lengthDevice, 1 * sizeof(int)); // keep track of the length of e_aff;
//    cudaMalloc((void**) &edgeCountDevice, 1 * sizeof(int)); // keep track of the length of E after long update;

    cudaMemcpy(edgeSrcDevice, edgeSrc, numEdges * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(edgeDstDevice, edgeDst, numEdges * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(rowPtrDevice,  rowPtr,  numRow * sizeof(int),   cudaMemcpyHostToDevice);


    dim3 dimBlock1(BLOCK_SIZE, 1, 1);
    dim3 dimGrid1(ceil(1.0*numEdges/BLOCK_SIZE), 1, 1);

//    flagInit<<<dimGrid1, dimBlock1>>>(edgeFlagDevice, numEdges);
    cudaDeviceSynchronize();

    while (true) {
        dim3 dimBlock(BLOCK_SIZE, 1, 1);
        dim3 dimGrid(ceil(1.0*numEdges/BLOCK_SIZE), 1, 1);

        // first mark all edge available to be affected and un-deleted
        mark<<<dimGrid, dimBlock>>>(edgeSrcDevice, affectedDevice, e_affDevice, to_deleteDevice, numEdges);
        cudaDeviceSynchronize();
//        printf("Pre-processing all data!\n");

        while (true) {

            // select e_aff out of affected
            length.data()[0] = 0;
            selectAff<<<dimGrid, dimBlock>>>(e_affDevice, affectedDevice, numEdges, length.data());
            cudaDeviceSynchronize();
//            printf("select aff!\n");

            if (length.data()[0] == 0) {
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
//            shortUpdate<<<dimGrid, dimBlock>>>(edgeSrcDevice, edgeDstDevice, to_deleteDevice, edgeFlagDevice, numEdges);
            shortUpdate<<<dimGrid, dimBlock>>>(edgeSrcDevice, edgeDstDevice, to_deleteDevice, numEdges);
            cudaDeviceSynchronize();

//            printf("Short update done!\n");
        }

        // long update
//        int * SDevice;
//
//        cudaMalloc((void**) &SDevice, ceil(numEdges/(2.0*BLOCK_SIZE)) * sizeof(int));
////
//        dim3 dimGrid2(ceil(numEdges/(2.0*BLOCK_SIZE)), 1, 1);
//        scan<<<dimGrid2, dimBlock>>>(edgeFlagDevice, edgePreSumDevice, SDevice, numEdges);
//        cudaDeviceSynchronize();
////
//        int numBlocks = ceil(numEdges/(2.0*BLOCK_SIZE));
//        dim3 dimGrid3(ceil(numBlocks/(2.0*BLOCK_SIZE)), 1, 1);
//        scan2<<<dimGrid3, dimBlock>>>(SDevice, SDevice, numBlocks);
//        cudaDeviceSynchronize();
////
//        dim3 dimBlock2(BLOCK_SIZE*2, 1, 1);
//        addScanBlockSum<<<dimGrid2, dimBlock2>>>(SDevice, edgePreSumDevice, numEdges);
//        cudaDeviceSynchronize();

        edgeCount.data()[0] = 0;
        countEdges<<<dimGrid, dimBlock>>>(edgeSrcDevice, numEdges, edgeCount.data());
        cudaDeviceSynchronize();

//        printf("Long update done!\n");

        if (edgeCount.data()[0] > 0) {
            printf("Order of truss is %d. Remaining # of edges is %d\n", k, edgeCount.data()[0]);
            k += 1;
        } else {
            k -= 1;
            break;
        }

//        cudaFree(SDevice);
    }

    cudaFree(edgeSrcDevice);
    cudaFree(edgeDstDevice);
    cudaFree(rowPtrDevice);
    cudaFree(affectedDevice);
    cudaFree(to_deleteDevice);
    cudaFree(e_affDevice);
//    cudaFree(edgeFlagDevice);
//    cudaFree(edgePreSumDevice);

//    cudaFree(lengthDevice);
//    cudaFree(edgeCountDevice);

    return k;
}

int do_k_truss(pangolin::COOView<int> view) {

//    // test case for our own
//    printf("For the first mini-test graph:\n");
//    int numEdges = 16;
//    int edgeSrc[16] = {0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};
//    int edgeDst[16] = {1, 2, 0, 2, 3, 4, 0, 1, 3, 4, 1, 2, 4, 1, 2, 3};
//    int rowPtr[6] = {0, 2, 6, 10, 13, 16};
//    int numRow = 6;
//    printf("k-truss order for the first graph is: %d\n", k_truss(edgeSrc, edgeDst, rowPtr, numEdges, numRow));
//
//    printf("For the second mini-test graph:\n");
//    int numEdges2 = 30;
//    int edgeSrc2[30] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8};
//    int edgeDst2[30] = {1, 2, 7, 0, 3, 4, 0, 5, 6, 7, 1, 4, 8, 1, 3, 5, 8, 2, 4, 6, 2, 5, 7, 8, 0, 2, 6, 3, 4, 6};
//    int rowPtr2[10] = {0, 3, 6, 10, 13, 17, 20, 24, 27, 30};
//    int numRow2 = 10;
//    printf("k-truss order for the second graph is %d\n", k_truss(edgeSrc2, edgeDst2, rowPtr2, numEdges2, numRow2));
//
    printf("For the real input graph - California Road Network:\n");
    int numEdges3 = view.nnz();
    int * edgeSrc3 = const_cast<int*>(view.row_ind());
    int * edgeDst3 = const_cast<int*>(view.col_ind());
    int * rowPtr3 = const_cast<int*>(view.row_ptr());
    int numRow3 = view.num_rows();
    printf("number of edges: %d\nnumber of nodes: %d\n", numEdges3, numRow3);
    printf("k-truss order for the input graph is %d\n", k_truss(edgeSrc3, edgeDst3, rowPtr3, numEdges3, numRow3));
    return 0;
}
