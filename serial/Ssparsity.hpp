#include "../include/utility.h"
#include<iostream>
#include<algorithm>
#include<string>
#include<math.h>
#include<stdlib.h>
#include<string.h>


class Ssparsity
{
public:
    int *deg = NULL;
    double *vol = NULL;
    int *ai = NULL, *aj = NULL, *insparsifier = NULL, SIZE, M, maxnode, mycnt;
    int *layer = NULL;
    double *av = NULL;
    DisjSet dis;
    Graph gra;
    bfsnode *bfsque;
    Node *mynum = NULL;

    void RadixSort(Node *num, int len)
    {
        Node *b = (Node *)malloc(sizeof(Node) * len);
        int sum[8][256]; memset(sum, 0, 8 * 256 * sizeof(int));
        
        for(int i = 0; i < len; i++)
        {
            for(int j = 0; j < 8; j++)
                ++sum[j][(reinterpret_cast<unsigned long long int&>(num[i].value) >> (j * 8)) & 255];
        }
        for(int q = 1; q <= 255; q++)
        {
            for(int j = 0; j < 8; j++) sum[j][q] += sum[j][q - 1];
        }
        for(int i = 0; i < 8; i++)
        {
            std::swap(num, b);
            int q = len - 1;
            while(1)
            {
                num[--sum[i][(reinterpret_cast<unsigned long long int&>(b[q].value) >> (i * 8)) & 255]] = b[q];
                if(q == 0) break;
                q--;
            }
        }
        free(b);
    }
    void initialize()
    {
        mycnt = 0;
        dis.func(SIZE);
        mynum = (Node *)malloc(sizeof(Node) * M);
        layer = (int*)malloc(sizeof(int) * (SIZE));
        bfsque = (bfsnode*)malloc(sizeof(bfsnode) * (SIZE));
        deg = (int*)malloc(sizeof(int) * (SIZE));
        vol = (double*)malloc(sizeof(double) * (SIZE));
        memset(insparsifier, 0, sizeof(int) * M);
        memset(layer, -1, sizeof(int) * SIZE);
        memset(deg, 0, sizeof(int) * SIZE);
        memset(vol, 0, sizeof(double) * SIZE);
        for (int l = 0; l < M; l++) mynum[l].idx = l;
    }
    void constructgraph()
    {
        gra = Graph(SIZE, M);
        gra.ptr = (int *)malloc(sizeof(int) * (SIZE + 1));
        gra.edgelists = (int *)malloc(sizeof(int) * (2 * M));
        int *ptrs = (int *)malloc(sizeof(int) * (SIZE));
        int i, j;
        for (int l = 0; l < M; l++)
        {
            i = ai[l]; j = aj[l];
            ++deg[i];
            ++deg[j];
            vol[i] += av[l];
            vol[j] += av[l];
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
        for (int l = 0; l < M; l++)
        {
            i = ai[l]; j = aj[l];
            gra.edgelists[ptrs[i]] = j;
            gra.edgelists[ptrs[j]] = i;
            ++ptrs[i]; ++ptrs[j];
        }
        free(ptrs);
    }
    void BFS()
    {
        int head = 0, tail = 0, node, curlayer, i = maxnode, tar;
        bfsque[tail] = bfsnode{ i,0 };
        tail++;
        layer[i] = 0;
        while (head < tail)
        {
            node = bfsque[head].node;
            curlayer = bfsque[head].layer;
            ++head;
            for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; ++k)
            {
                tar = gra.edgelists[k];
                if (layer[tar] < 0)
                {
                    bfsque[tail++] = bfsnode{ tar, curlayer + 1 };
                    layer[tar] = curlayer + 1;
                }
            }
        }
        for(int i = 0; i < SIZE; i++)
        {
            if(layer[i] >= 0) continue;
            head=0;tail=0;
            bfsque[tail]=(bfsnode){i,0};tail++;
            layer[i]=0;
            while(head < tail)
            {
                node=bfsque[head].node; curlayer=bfsque[head].layer; head++;
                for (int k = gra.ptr[node]; k < gra.ptr[node + 1]; ++k)
                {
                    tar = gra.edgelists[k];
                    if (layer[tar] < 0)
                    {
                        bfsque[tail++] = bfsnode{ tar, curlayer + 1 };
                        layer[tar] = curlayer + 1;
                    }
                }
            }
        }
    }
    void Kruskal()
    {
        double maxdeg, scale = 1;
        int i, j, l, ind, count = 0;
        for (l = 0; l < M; l++)
        {
            i = ai[l]; j = aj[l];
            maxdeg = std::max(deg[i], deg[j]);
            scale = (maxdeg) / double(log((layer[i] + layer[j]) + 1));
            mynum[l].value = av[l] * scale; 
        }
        RadixSort(mynum, M);
        l = M-1;
        while(1)
        {
            ind = mynum[l].idx;
            i = ai[ind]; j = aj[ind];
            if (!dis.is_same(i, j))
            {
                insparsifier[ind] = true;
                ++count;
                if(count >= SIZE - 1) break;
                dis.to_union(i, j);
            }
            if(l <= 0) break;
            --l;
        }
        for(int a = 0; a < M; ++a)
        {
            if(!insparsifier[a])
            {
                mynum[mycnt].idx = a;
                mycnt++;
            }
        }
    }
    void calculateresistance()
    {
        int i, j, k, p, l;
        for (int a = 0; a < mycnt; a++)
        {
            l = mynum[a].idx;
            j = ai[l];
            k = aj[l];
            mynum[a].value = av[l] / (vol[j] + vol[k] + 2 * av[l]);
        }
        RadixSort(mynum, mycnt);
    }
    void addedges(double fac)
    {
        int num = int(fac * double(SIZE)), cnt = 0; if(num > mycnt) num = mycnt;
        for(int l = mycnt - num, i; l < mycnt; ++l)
        {
            i = mynum[l].idx;
            if(insparsifier[i]) ++cnt;
            insparsifier[i] = true;
        }
    }
    void freeMemory()
    {
        dis.destroy();
        gra.destroy();
        free(deg); free(layer); free(bfsque); free(vol);
        free(mynum);
    }
    void Sparsity(int *ai_in, int *aj_in, double *av_in,
        int M_in, int N_in, int *insparsifier_in, double alpha_in)
    {
        //struct timeval t1, t2;
        ai = ai_in; aj = aj_in; av = av_in; M = M_in; SIZE = N_in; insparsifier = insparsifier_in; double fac = alpha_in;
        initialize();
        constructgraph();
        //gettimeofday(&t1, NULL);
        BFS();
        Kruskal();
        calculateresistance();
        addedges(fac);
        //gettimeofday(&t2, NULL);
        freeMemory();
    }
    void Sparsify(double p, int n, int m, int &m2,
        int *i1, int *j1, double *w1,
        int **i2, int **j2, double **w2)
    {
        int *i = (int *)malloc(sizeof(int) * m);
        int *j = (int *)malloc(sizeof(int) * m);
        int *in = (int *)malloc(sizeof(int) * m); memset(in, 0, sizeof(int) * m);
        double *w = (double *)malloc(sizeof(double) * m);
        for (int a = 0; a < m; a++)
        {
            i[a] = i1[a] - 1;
            j[a] = j1[a] - 1;
            w[a] = fabs(w1[a]);
        }
        Sparsity(i, j, w, m, n, in, p);
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
        free(i);
        free(j);
        free(w);
        free(in);
        return;
    }
};
