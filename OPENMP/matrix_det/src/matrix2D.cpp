#include "../include/matrix2D.h"
#include <stdlib.h>

matrix2D::matrix2D(int n, int m)
{
    this->mat = (int *)calloc(n * m, sizeof(int));
    this->n = n;
    this->m = m;
}

void matrix2D::set(int x, int y, int value)
{
    int i = this->n * x + y;
    this->mat[i] = value;
}

int matrix2D::get(int x, int y)
{
    if (x >= n || y >= m)
        return 0;

    int i = this->n * x + y;
    return this->mat[i];
}

int matrix2D::get_n()
{
    return this->n;
}

int matrix2D::get_m()
{
    return this->m;
}