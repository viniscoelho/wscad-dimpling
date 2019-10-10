#include <cmath>
#include <cstdio>
#include <cstdlib>

#define C 4
#define THREADS 1024 // 2^10
#define MAX 85
#define MAXS MAX* MAX
#define COMB_MAX (MAX * (MAX - 1) * (MAX - 2) * (MAX - 3)) / 24

#define gpuErrChk(ans)                        \
    {                                         \
        gpuAssert((ans), __FILE__, __LINE__); \
    }
inline void gpuAssert(cudaError_t code, char* file, int line, bool abort = true)
{
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),
            file, line);
        if (abort)
            getchar();
    }
}

using namespace std;

struct Node {
    int sz, qtd;
    int graph[MAXS], TMP[6 * MAX], seeds[COMB_MAX * C];
};

struct Params {
    int faces, count, tmpMax;
    int F[6 * MAX], V[MAX];
};

/*
    SIZE        ---> Number of vertices
    C           ---> Size of the combination (Size of a seed clique)
    faces       ---> Quantity of triangular faces
    qtd         ---> Number of possible 4-cliques
    T           ---> Output graph for an instance
    R           ---> Output graph for an possible optimal solution
    F           ---> List containing triangular faces of an instance
    seeds       ---> Combinations of possible starting 4-cliques
    graph       ---> The graph itself
*/
double start, stop;
int R[MAX][MAX], F[8 * MAX], bib[MAX];
int SIZE, BLOCKS, COMB, qtd = 0;

Node* N;

//-----------------------------------------------------------------------------
// Mac OSX
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

/*
    Prints elapsed time.
    */
void printElapsedTime(double start, double stop)
{
    double elapsed = stop - start;
    printf("Elapsed time: %.3lfs.\n", elapsed);
}
//-----------------------------------------------------------------------------
/*  
    Gets clock time.
    */
void current_utc_time(struct timespec* ts)
{
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    ts->tv_sec = mts.tv_sec;
    ts->tv_nsec = mts.tv_nsec;
#else
    clock_gettime(CLOCK_REALTIME, ts);
#endif
}
//-----------------------------------------------------------------------------
double getTime()
{
    timespec ts;
    current_utc_time(&ts);
    return double(ts.tv_sec) + double(ts.tv_nsec) / 1e9;
}
//-----------------------------------------------------------------------------
/*
    t   ---> thread index
    Generates a list of vertices which are not on the initial planar graph.
    */
__device__ void generateVertexList(Node* devN, Params* devP, int t)
{
    int sz = devN->sz;
    int va = devN->seeds[t * 4], vb = devN->seeds[t * 4 + 1],
        vc = devN->seeds[t * 4 + 2], vd = devN->seeds[t * 4 + 3];
    for (int i = 0; i < sz; ++i) {
        if (i == va || i == vb || i == vc || i == vd)
            devP[t].V[i] = -1;
        else
            devP[t].V[i] = i;
    }
}
//-----------------------------------------------------------------------------
/*
    t   ---> thread index
    Returns the initial solution weight for the planar graph and
    initializes necessary structures, such as the edges indexes,
    and defines which edges belongs to a face.
    */
__device__ void generateFaceList(Node* devN, Params* devP, int t)
{
    int resp = 0, sz = devN->sz;
    int va = devN->seeds[t * 4], vb = devN->seeds[t * 4 + 1],
        vc = devN->seeds[t * 4 + 2], vd = devN->seeds[t * 4 + 3];

    // Generate the first triangle of the output graph
    devP[t].F[devP[t].faces * 3] = va,
                              devP[t].F[devP[t].faces * 3 + 1] = vb,
                              devP[t].F[(devP[t].faces++) * 3 + 2] = vc;
    resp = devN->graph[va * sz + vb] + devN->graph[va * sz + vc] + devN->graph[vb * sz + vc];

    // Generate the next 3 possible faces
    devP[t].F[devP[t].faces * 3] = va, devP[t].F[devP[t].faces * 3 + 1] = vb,
                              devP[t].F[(devP[t].faces++) * 3 + 2] = vd;
    devP[t].F[devP[t].faces * 3] = va, devP[t].F[devP[t].faces * 3 + 1] = vc,
                              devP[t].F[(devP[t].faces++) * 3 + 2] = vd;
    devP[t].F[devP[t].faces * 3] = vb, devP[t].F[devP[t].faces * 3 + 1] = vc,
                              devP[t].F[(devP[t].faces++) * 3 + 2] = vd;
    resp += devN->graph[va * sz + vd] + devN->graph[vb * sz + vd] + devN->graph[vc * sz + vd];
    devP[t].tmpMax = resp;
}
//-----------------------------------------------------------------------------
/*
    Inserts a new vertex, 3 new triangular faces
    and removes the face from the list.
    */
__device__ int faceDimple(Node* devN, Params* devP, int new_vertex, int f, int t)
{
    // Remove the chosen face and insert a new one
    int va = devP[t].F[f * 3],
        vb = devP[t].F[f * 3 + 1],
        vc = devP[t].F[f * 3 + 2];

    devP[t].F[f * 3] = new_vertex,
                  devP[t].F[f * 3 + 1] = va, devP[t].F[f * 3 + 2] = vb;
    // Insert the other two possible faces
    devP[t].F[devP[t].faces * 3] = new_vertex,
                              devP[t].F[devP[t].faces * 3 + 1] = va,
                              devP[t].F[(devP[t].faces++) * 3 + 2] = vc;
    devP[t].F[devP[t].faces * 3] = new_vertex,
                              devP[t].F[devP[t].faces * 3 + 1] = vb,
                              devP[t].F[(devP[t].faces++) * 3 + 2] = vc;

    int sz = devN->sz;
    int resp = devN->graph[va * sz + new_vertex] + devN->graph[vb * sz + new_vertex]
        + devN->graph[vc * sz + new_vertex];

    return resp;
}
//-----------------------------------------------------------------------------
/*
    Returns the vertex having the maximum gain
    inserting within a face.
    */
__device__ int maxGainFace(Node* devN, Params* devP, int* f, int t)
{
    int sz = devN->sz;
    int gain = -1, vertex = -1;

    int faces = devP[t].faces;
    // Iterate through the remaining vertices
    for (int new_vertex = 0; new_vertex < sz; ++new_vertex) {
        if (devP[t].V[new_vertex] == -1)
            continue;
        // Test the dimple on each face
        for (int i = 0; i < faces; ++i) {
            int va = devP[t].F[i * 3], vb = devP[t].F[i * 3 + 1], vc = devP[t].F[i * 3 + 2];
            int tmpGain = devN->graph[va * sz + new_vertex] + devN->graph[vb * sz + new_vertex]
                + devN->graph[vc * sz + new_vertex];
            if (tmpGain > gain) {
                gain = tmpGain;
                *f = i;
                vertex = new_vertex;
            }
        }
    }
    return vertex;
}
//-----------------------------------------------------------------------------
__device__ void dimpling(Node* devN, Params* devP, int t)
{
    while (devP[t].count) {
        int f = -1;
        int vertex = maxGainFace(devN, devP, &f, t);
        devP[t].V[vertex] = -1;
        devP[t].tmpMax += faceDimple(devN, devP, vertex, f, t);
        devP[t].count--;
    }
}
//-----------------------------------------------------------------------------
__device__ void initializeDevice(Params* devP, int sz, int t)
{
    devP[t].faces = 0;
    devP[t].tmpMax = -1;
    devP[t].count = sz - 4;
}
//-----------------------------------------------------------------------------
__global__ void solve(Node* devN, Params* devP, int* respMax, int* idx)
{
    int x = blockDim.x * blockIdx.x + threadIdx.x;
    int sz = devN->sz;
    int comb = devN->qtd;
    __syncthreads();

    if (x < comb) {
        initializeDevice(devP, devN->sz, x);
        generateVertexList(devN, devP, x);
        generateFaceList(devN, devP, x);
        dimpling(devN, devP, x);

        __syncthreads();
        atomicMax(respMax, devP[x].tmpMax);
        if (devP[x].tmpMax == *respMax)
            *idx = x;
        __syncthreads();
    }
}
//-----------------------------------------------------------------------------
int prepare()
{
    int resp = 0, idx = 0, *tmpResp, *tmpIdx;
    gpuErrChk(cudaMalloc((void**)&tmpResp, sizeof(int)));
    gpuErrChk(cudaMalloc((void**)&tmpIdx, sizeof(int)));
    gpuErrChk(cudaMemcpy(tmpResp, &resp, sizeof(int), cudaMemcpyHostToDevice));
    gpuErrChk(cudaMemcpy(tmpIdx, &idx, sizeof(int), cudaMemcpyHostToDevice));

    Node* devN;
    Params* devP;

    size_t sz = COMB * sizeof(Params);
    gpuErrChk(cudaMalloc((void**)&devP, sz));
    gpuErrChk(cudaMalloc((void**)&devN, sizeof(Node)));
    gpuErrChk(cudaMemcpy(devN, N, sizeof(Node), cudaMemcpyHostToDevice));

    dim3 blocks(BLOCKS, 1);
    dim3 threads(THREADS, 1);

    solve<<<blocks, threads>>>(devN, devP, tmpResp, tmpIdx);
    gpuErrChk(cudaDeviceSynchronize());

    gpuErrChk(cudaMemcpy(&resp, tmpResp, sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrChk(cudaMemcpy(&idx, tmpIdx, sizeof(int), cudaMemcpyDeviceToHost));
    gpuErrChk(cudaMemcpy(&F, devP[idx].F, (6 * MAX) * sizeof(int), cudaMemcpyDeviceToHost));

    gpuErrChk(cudaFree(devN));
    gpuErrChk(cudaFree(devP));
    return resp;
}
//-----------------------------------------------------------------------------
/*
    C      ---> Size of the combination
    index  ---> Current index in data[]
    data[] ---> Temporary array to store a current combination
    i      ---> Index of current element in vertices[]
*/
void combineUntil(int index, int* data, int i)
{
    if (index == C) {
        for (int j = 0; j < C; ++j)
            N->seeds[qtd * C + j] = data[j];
        qtd++;
        return;
    }

    if (i >= SIZE)
        return;

    data[index] = i;
    combineUntil(index + 1, data, i + 1);
    combineUntil(index, data, i + 1);
}
//-----------------------------------------------------------------------------
void combine()
{
    int data[C];
    combineUntil(0, data, 0);
}
//-----------------------------------------------------------------------------
/*
    Defines the number of combinations.
    */
void sizeDefinitions()
{
    for (int i = 4; i <= MAX; ++i) {
        int resp = 1;
        for (int j = i - 3; j <= i; ++j)
            resp *= j;
        resp /= 24;
        bib[i - 1] = resp;
    }
}
//-----------------------------------------------------------------------------
void initialize()
{
    for (int i = 0; i < SIZE; ++i) {
        for (int j = i + 1; j < SIZE; ++j)
            R[i][j] = R[j][i] = -1;
        R[i][i] = -1;
    }
}
//-----------------------------------------------------------------------------
void readInput()
{
    int x;
    scanf("%d", &SIZE);
    COMB = bib[SIZE - 1];
    BLOCKS = COMB / THREADS + 1;

    N = (Node*)malloc(sizeof(Node));
    N->sz = SIZE;
    N->qtd = COMB;

    for (int i = 0; i < SIZE; ++i) {
        for (int j = i + 1; j < SIZE; ++j) {
            scanf("%d", &x);
            N->graph[i * SIZE + j] = x;
            N->graph[j * SIZE + i] = x;
        }
        N->graph[i * SIZE + i] = -1;
    }
}
//-----------------------------------------------------------------------------
int main(int argv, char** argc)
{
    sizeDefinitions();
    // Reads the input, which is given by the size of a graph and its weighted
    // edges. The given graph should be a complete graph.
    readInput();
    initialize();
    // Generate 4-clique seeds, given the number of vertices
    combine();

    start = getTime();
    int respMax = prepare();
    stop = getTime();

    for (int i = 0; i < 2 * SIZE; ++i) {
        int va = F[i * 3], vb = F[i * 3 + 1], vc = F[i * 3 + 2];
        if (va == vb && vb == vc)
            continue;
        R[va][vb] = R[vb][va] = N->graph[va * SIZE + vb];
        R[va][vc] = R[vc][va] = N->graph[va * SIZE + vc];
        R[vb][vc] = R[vc][vb] = N->graph[vb * SIZE + vc];
    }

    printf("Printing generated graph:\n");
    for (int i = 0; i < SIZE; ++i) {
        for (int j = i + 1; j < SIZE; ++j)
            printf("%d ", R[i][j]);
        printf("\n");
    }

    printElapsedTime(start, stop);
    printf("Maximum weight found: %d\n", respMax);
    free(N);

    return 0;
}
