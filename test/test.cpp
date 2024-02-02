#include "../serial/Ssparsity.hpp"
#include "../parallel/Psparsity.hpp"
#include <cstdio>
#include <cstring>
#include <cmath>
void ReadData(char *path, int &n1, int &m1, int **i1, int **j1, double **w1) // 1.文件绝对路径 3.边i 4.边j 5.权重w
{
    FILE *fp; fp = fopen(path, "r");
    unsigned temp1; int temp2;
    temp2 = fscanf(fp, "%u %u %u", &n1, &temp1, &m1);
    *i1 = (int *)malloc(sizeof(int) * m1);
    *j1 = (int *)malloc(sizeof(int) * m1);
    *w1 = (double *)malloc(sizeof(double) * m1);
    int ii, jj, cnt = 0, cnt1 = 0, mm = 0;
    double ww;
    for(cnt = 0; cnt < m1; ++cnt)
    {
        temp2 = fscanf(fp, "%u %u %lf", &ii, &jj, &ww);
        if(jj < ii)
        {
            (*i1)[cnt1] = ii;
            (*j1)[cnt1] = jj;
            (*w1)[cnt1] = ww;
            ++cnt1;
        }
    }
    fclose(fp);
    m1 = cnt1;
}
int main(int argc, char **argv)
{
    int *i1 = NULL, *j1 = NULL, *i2 = NULL, *j2 = NULL, *i3 = NULL, *j3 = NULL, n, m, m2, m3;
    double *w1, *w2 = NULL, *w3 = NULL;
    ReadData(argv[1], n, m, &i1, &j1, &w1);
    Ssparsity sparsity1;
    sparsity1.Sparsify(atof(argv[2]), n, m, m2, i1, j1, w1, &i2, &j2, &w2);
    printf("n = %-10d    edges1 = %-10d    edges2 = %-10d\n", n, m, m2);
    free(i2); i2 = NULL;
    free(j2); j2 = NULL;
    free(w2); w2 = NULL;
    //(n, m, i1, j1, w1, 0.02, 8, m2, &i2, &j2, &w2);
    Psparsity sparsity2;
    sparsity2.Sparsify(atof(argv[2]), n, m, m2, atoi(argv[3]), i1, j1, w1, &i2, &j2, &w2);
    printf("n = %-10d    edges1 = %-10d    edges2 = %-10d\n", n, m, m2);
    free(i1);
    free(j1);
    free(w1);
    free(j2);
    free(i2);
    free(w2);
    return 0;
}