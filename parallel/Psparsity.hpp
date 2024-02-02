#include "../include/utility.h"

#include <atomic>
#include <iostream>
#include <algorithm>
#include <execution>
#include <chrono>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <omp.h>


class PBfs
{
public:
    int SIZE, nthreads, globalLen, level;
    int *globalQueue, *prefixSum, **localQueuesList;
    int *layer, *ptr;
    CSREdge *edgelists;
    PBfs(int SIZE, int nthreads, int *layer, int *ptr, CSREdge *edgelists):
    globalQueue(NULL), prefixSum(NULL), localQueuesList(NULL), SIZE(SIZE), nthreads(nthreads), layer(layer), ptr(ptr), edgelists(edgelists)
    {
        globalQueue = (int *)malloc(sizeof(int) * SIZE);
        prefixSum = (int *)malloc(sizeof(int) * (nthreads + 1));
        localQueuesList = (int **)malloc(sizeof(int *) * (nthreads));
        for(int a = 0; a < nthreads; ++a) localQueuesList[a] = (int *)malloc(sizeof(int) * SIZE);
    }
    bool PBFS1()
    {
        bool improvement = false;

        #pragma omp parallel num_threads(nthreads)
        {
            int tid = omp_get_thread_num();
            int *localQueue = localQueuesList[tid];
            int localLen = 0, p1, s, e, u;

            #pragma omp for reduction(||:improvement) schedule(guided, 32)
            for (int i = 0; i < globalLen; i++)
            {
                p1 = globalQueue[i];
                s = ptr[p1]; e = ptr[p1 + 1];
                for (; s < e; s++)
                {
                    u = edgelists[s].node;
                    if(layer[u] < 0)
                    {
                        layer[u] = level + 1;
                        localQueue[localLen++] = u;
                        improvement = true;
                    }
                }
            }
            
            prefixSum[tid + 1] = localLen;
            #pragma omp barrier
            #pragma omp single
            {
                for (int i = 0; i < omp_get_num_threads(); i++) prefixSum[i + 1] += prefixSum[i];
            }

            memcpy(globalQueue + prefixSum[tid], localQueue, sizeof(int) * (prefixSum[tid + 1] - prefixSum[tid]));
            globalLen = prefixSum[omp_get_num_threads()]; 
        }

        if(improvement) ++level;
        return improvement;
    }
    void PBFS(int root)
    {
        layer[root] = 0;
        globalQueue[0] = root;
        prefixSum[0] = 0;
        globalLen = 1, level = 0;
        while(PBFS1()) ;
    }
    void destroy()
    {
        free(globalQueue);
        free(prefixSum);
        for(int a = 0; a < nthreads; ++a) free(localQueuesList[a]);
    }
};
class Psparsity
{
public:
    int *deg = NULL, *ptrs = NULL, *insparsifier = NULL, *ptr1 = NULL;
    double *vol = NULL;
    int SIZE, M, maxnode, len1, len2, mycnt, NTHREADS;
    int *layer = NULL, *intNum1;
    std::atomic_uint atomicuint1 = 0;
    Edge *edges = NULL;
    DisjSet dis;
    PSGraph gra;
    Node *mynum = NULL;

    void initialize()
    {
        len1 = SIZE / NTHREADS;
        len2 = M / NTHREADS;
        mycnt = 0;

        gra = PSGraph(SIZE, M);
        gra.ptr = (int *)malloc(sizeof(int) * (SIZE + 1));
        gra.edgelists = (CSREdge *)malloc(sizeof(CSREdge) * (2 * M));
        intNum1 = (int*)malloc(sizeof(int) * (SIZE)); memset(intNum1, -1, sizeof(int) * (SIZE));
        ptrs = (int *)malloc(sizeof(int) * (SIZE));
        ptr1 = (int *)malloc(sizeof(int) * (SIZE+1));
        dis.func(SIZE);
        mynum = (Node *)malloc(sizeof(Node) * M);
        layer = (int*)malloc(sizeof(int) * (SIZE));
        deg = (int *)malloc(sizeof(int) * SIZE);
        vol = (double *)malloc(sizeof(double) * SIZE);

        memset(dis.rank, 0, sizeof(int) * SIZE);
        for(int i = 0; i < SIZE; ++i) dis.parent[i] = i;
        memset(insparsifier, 0, sizeof(int) * M);
        memset(layer, -1, sizeof(int) * SIZE);
        for(int i = 0; i < SIZE; ++i) deg[i] = 0;
        for(int i = 0; i < SIZE; ++i) vol[i] = 0;
        for(int i = 0; i < M; ++i) mynum[i].idx = i;
    }
    void constructgraph()
    {
        int i, j, temp1, p;
        double w;
        for(int l = 0; l < M; ++l)
        {
            i = edges[l].u; j = edges[l].v; w = edges[l].w;
            ++deg[i];
            ++deg[j];
            vol[i] += w;
            vol[j] += w;
        }
        maxnode = 0; gra.ptr[0] = 0; ptrs[0] = 0;
        int sum = 0;
        for (int l = 1; l < SIZE; l++)
        {
            if (vol[l] > vol[maxnode]) maxnode = l;
            sum += deg[l - 1];
            gra.ptr[l] = sum;
            ptrs[l] = sum;
        }
        sum += deg[SIZE - 1];
        gra.ptr[SIZE] = sum;
        for (int l = 0; l < M; ++l)
        {
            i = edges[l].u; j = edges[l].v;
            gra.edgelists[ptrs[i]] = CSREdge{j, l};
            gra.edgelists[ptrs[j]] = CSREdge{i, l};
            ++ptrs[i]; ++ptrs[j];
        }
    }
    void caleffw(int &tid)
    {
        int s = len2 * tid, e, i, j, l;
        double maxdeg;
        if(tid == NTHREADS - 1) e = M;
        else e = s + len2;
        for(; s < e; ++s)
        {
            i = edges[s].u; j = edges[s].v;
            maxdeg = std::max(deg[i], deg[j]);
            mynum[s].value = edges[s].w * ((maxdeg) / double(log((layer[i] + layer[j]) + 1)));
        }
    }
    void ST1(int &tid)
    {
        int s, e;
        s = len1 * tid;
        if(tid == NTHREADS - 1) e = SIZE;
        else e = s + len1;
        for(; s < e; ++s)
            if(gra.ptr[s+1] - gra.ptr[s] == 1) insparsifier[gra.ptr[s]] = true;
    }
    void ST2(int &tid)
    {
        CSREdge *pedges = gra.edgelists, *pedge, *end;
        int tar, tar1 = len1 * tid, *pptr = gra.ptr, s, e, r1, end1;
        double pmax, val;
        if(tid == NTHREADS - 1) end1 = SIZE;
        else end1 = tar1 + len1;
        for(; tar1 < end1; ++tar1)
        {
            s = pptr[tar1], e = pptr[tar1 + 1];
            pedge = &pedges[s], end = &pedges[e];
            tar = r1 = pedge->ind;
            if(e - s <= 1) continue;
            pmax = mynum[tar].value, pedge++;
            for(; pedge != end; pedge++)
            {
                tar = pedge->ind;
                val = mynum[tar].value;
                if(pmax <= val)
                {
                    pmax = val;
                    r1 = tar;
                }
            }
            intNum1[tar1] = r1;
        }
    }
    void ST3(int &tid)
    {
        Edge *pEdge;
        int tar, tar1 = len1 * tid, end1, r1, r2;
        if(tid == NTHREADS - 1) end1 = SIZE;
        else end1 = tar1 + len1;
        for(; tar1 < end1; ++tar1)
        {
            tar = intNum1[tar1];
            if(tar == -1) continue;
            pEdge = &edges[tar];
            r1 = pEdge->u, r2 = pEdge->v;
            if(tar1 == r2 && intNum1[r1] == tar)
            {
                dis.to_union(r1, r2);
                *(insparsifier + tar) = true;
            }
        }
    }
    void kruskal()
    {
        int count = M - mycnt, i, j, ind;
        int l = mycnt-1;
        size_t size1 = sizeof(Node);
        while(1)
        {
            ind = mynum[l].idx;
            i = edges[ind].u; j = edges[ind].v;
            if(!dis.is_same(i, j))
            {
                insparsifier[ind] = true;
                ++count;
                if(count == SIZE - 1) break;
                dis.to_union(i, j);
            }
            if(l == 0) break;
            --l;
        }
        mycnt = 0;
        for(int a = 0; a < M; ++a)
        {
            if(!insparsifier[a])
            {
                mynum[mycnt].idx = a;
                mycnt++;
            }
        }
    }
    void calR(int &tid)
    {
        int len3 = mycnt / NTHREADS;
        int s = len3 * tid, e, j, k, l;
        if(tid == NTHREADS - 1) e = mycnt;
        else e = s + len3;
        for(; s < e; ++s)
        {
            l = mynum[s].idx;
            j = edges[l].u;
            k = edges[l].v;
            mynum[s].value = edges[l].w / (vol[j] + vol[k] + 2 * edges[l].w);
        }
    }
    void addedges(double &fac, int &tid)
    {
        int i, num = int(fac * double(SIZE)), s, e;
        s = mycnt - num / NTHREADS * tid - 1;
        if(tid == NTHREADS - 1) e = mycnt - num;
        else e = mycnt - num / NTHREADS * (tid + 1) - 2;
        for(; s >= e; --s)
        {
            i = mynum[s].idx;
            insparsifier[i] = true;
        }
    }
    void freeMemory()
    {
        dis.destroy();
        gra.destroy();
        free(deg); free(layer); free(vol); free(intNum1);
        free(mynum);
        free(ptrs);
        free(ptr1);
    }
    void sparsity(Edge *edges_in, int M_in, int N_in, int nthreads
        , int *insparsifier_in, double alpha_in)
    {
        //struct timeval t1, t2, t3, t4;
        edges = edges_in; M = M_in; SIZE = N_in; insparsifier = insparsifier_in;
        double fac = alpha_in;
        NTHREADS = nthreads;
        omp_set_num_threads(NTHREADS);
        
        initialize();
        constructgraph();
        PBfs bfs(SIZE, NTHREADS, layer, gra.ptr, gra.edgelists); 
        //gettimeofday(&t1, NULL);
        bfs.PBFS(maxnode);
        for(int a = 0; a < SIZE; ++a)
        {
            if(layer[a] >= 0 || deg[a] == 0) continue;
            bfs.PBFS(a);
        }
        #pragma omp parallel num_threads(NTHREADS)
        {
            int tid = omp_get_thread_num();
            caleffw(tid);
            #pragma omp barrier
            ST1(tid);
            #pragma omp barrier
            ST2(tid);
            #pragma omp barrier
            ST3(tid);
        }
        for(int a = 0; a < M; ++a)
        {
            if(!insparsifier[a])
            {
                mynum[mycnt].idx = a;
                mynum[mycnt].value = mynum[a].value;
                mycnt++;
            }
        }
        std::sort(std::execution::par, mynum, mynum + mycnt, [this](Node i, Node j){return i.value > j.value;});
        kruskal();
        #pragma omp parallel num_threads(NTHREADS)
        {
            int tid = omp_get_thread_num();
            calR(tid);
        }
        std::sort(std::execution::par, mynum, mynum + mycnt, [this](Node i, Node j){return i.value > j.value;});
        #pragma omp parallel num_threads(NTHREADS)
        {
            int tid = omp_get_thread_num();
            addedges(fac, tid);
        }
        //gettimeofday(&t2, NULL);
        bfs.destroy();
        freeMemory();
    }
    void Sparsify(double p, int n, int m, int &m2, int nthreads,
        int *i1, int *j1, double *w1,
        int **i2, int **j2, double **w2)
    {
        Edge *edges = (Edge *)malloc(sizeof(Edge) * m);
        int *in = (int *)malloc(sizeof(int) * m);  memset(in, 0, sizeof(int) * m);
        for (int i = 0; i < m; i++)
        {
            edges[i].u = i1[i] - 1;
            edges[i].v = j1[i] - 1;
            edges[i].w = fabs(w1[i]);
        }
        sparsity(edges, m, n, nthreads, in, p);
        m2 = 0;
        for(int a = 0; a < M; ++a)
        {
            if(in[a]) ++m2;
        }
        if(i2 != NULL)
        {
            (*i2) = (int *)malloc(sizeof(int) * m2);
            (*j2) = (int *)malloc(sizeof(int) * m2);
            (*w2) = (double *)malloc(sizeof(double) * m2);
            for(int a = 0, b = 0; a < M; ++a)
            {
                if(in[a])
                {
                    (*i2)[b] = i1[a];
                    (*j2)[b] = j1[a];
                    (*w2)[b] = w1[a];
                    ++b;
                }
            }
        }
        free(edges);
        free(in);
        return;
    }
};