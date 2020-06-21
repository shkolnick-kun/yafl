#include <stdio.h>
#include "cyexp_c.h"

void do_callback(cHelloSt * arg)
{
    if (arg->cy_wrapper)
    {
        arg->cy_wrapper(arg);
        return;
    }
    printf("\nFFFFFFUUUU~~~~!!!\n");
}

void hello(void)
{
    printf("\nHello cython!\n");
}

void print_mat(int nr, int nc, double * m)
{
    int i;
    int j;
    printf("\n");
    for (i = 0; i < nr; i++)
    {
        printf("[%f", m[nc*i]);
        for (j = 1; j < nc; j++)
        {
            printf(" %f", m[nc*i + j]);
        }
        printf("]\n");
    }
}

void zero_mat(int nr, int nc, double * m)
{
    int i;
    int j;
    for (i = 0; i < nr; i++)
    {
        for (j = 0; j < nc; j++)
        {
            m[nc*i + j] = 0.0;
        }
    }
}

void bye(void)
{
    printf("\nBye cython!\n");
}
