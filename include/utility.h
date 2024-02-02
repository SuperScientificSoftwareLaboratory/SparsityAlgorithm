#pragma once
#include <cstdlib>
#include<sys/time.h>

typedef struct edge {
    int u, v;
    double w;
}Edge;

typedef struct CSRedge
{
    int node, ind;
}CSREdge;

typedef struct node
{
    int idx;
    double value;
}Node;

typedef struct bfsnode
{
    int node, layer;
}bfsnode;

typedef struct dfsnode
{
    int node, from;
}dfsnode;
class Graph
{
public:
    int n, m, *ptr;
    int *edgelists;
    Graph(int n, int m) :n(n), m(m) { }
    Graph() { m = 0; }
    void destroy() { free(edgelists); free(ptr); }
};
class PSGraph
{
public:
    int n, m, *ptr;
    CSREdge *edgelists;
    PSGraph(int n, int m) :n(n), m(m) { }
    PSGraph() { m = 0; }
    void destroy() { free(edgelists); free(ptr); }
};
class DisjSet
{
public:
    int *parent, *rank;
    void func(int max_size)
    {
        parent = (int *)malloc(sizeof(int) * max_size);
        rank = (int *)malloc(sizeof(int) * max_size);
        for(int i = 0; i < max_size; i++) parent[i] = i;
    }
    int find(int &x)
    {
        int r = x;
        while(parent[r] != r) r = parent[r];
        int i = x, j;
        while(i != r)
        {
            j = parent[i];
            parent[i] = r;
            i = j;
        }
        return r;
    }
    void to_union(int &x1, int &x2)
    {
        x1 = find(x1);
        x2 = find(x2);
        if (x1 == x2) return;
        else if (rank[x1] < rank[x2]) parent[x1] = x2;
        else
        {
            parent[x2] = x1;
            if (rank[x1] == rank[x2]) ++rank[x1];
        }
    }
    bool is_same(int &e1, int &e2)
    {
        e1 = find(e1);
        e2 = find(e2);
        return e1 == e2;
    }
    void destroy() { free(parent); free(rank); }
};