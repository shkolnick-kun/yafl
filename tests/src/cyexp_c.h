#ifndef CYEXP_H
#define CYEXP_H

typedef struct _cHelloSt cHelloSt;

struct _cHelloSt {
    void (*cy_wrapper)(cHelloSt *);
    void * py_func;
    void * py_object;
};

void chello(void);

void print_mat(int nr, int nc, double * m);
void zero_mat(int nr, int nc, double * m);

#endif // CYEXP_H
