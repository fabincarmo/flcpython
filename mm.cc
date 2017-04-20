//#include "mm.hpp"
#include <locale.h>
#include <iomanip>
#include "mmio.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <tuple>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void bmap(Eigen::VectorXd v, char * filename)
{
    int nz = 0;
    int N = v.rows();
    int N_ = 1;
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    mm_set_symmetric(&matcode);
    FILE * f = fopen(filename, "w");
    mm_write_banner(f, matcode);
    for(int i=0; i<N; i++){
        for(int j=0; j<N_; j++){
            if ( v(i) == v(j) ){
                nz++;
            }
        }
        N_++;
    }
    N_ = 1;
    mm_write_mtx_crd_size(f, N, N, nz);
    for(int i=0; i<N; i++){
        for(int j=0; j<N_; j++){
            if ( v(i) == v(j) ){
                fprintf(f, "%d %d 1\n", i+1, j+1);
            }
        }
        N_++;
    }
    fclose(f);
    return;
}

void kmap(Eigen::VectorXd x, float var, float lim, char * filename)
{
    int N, i, j;
    int nz = 0;
    int N_ = 1;
    float v, aux;
    setlocale(LC_NUMERIC,"en_US");
    N = x.rows();
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    mm_set_symmetric(&matcode);
    FILE * f = fopen(filename, "w");
    mm_write_banner(f, matcode);

    for(i=0; i<N; i++){
        for(j=0; j<N_; j++){
            if ( exp(-pow(x(i)-x(j),2)/var) > lim ){
                nz++;
            }
        }
        N_++;
    }
    mm_write_mtx_crd_size(f, N, N, nz);
    N_ = 1;
    for(i=0; i<N; i++){
        for(j=0; j<N_; j++){
            if ( exp(-pow(x(i)-x(j),2)/var) > lim ){
                aux = exp(-pow(x(i)-x(j),2)/var);
                fprintf(f, "%d %d %10.3g\n", i+1, j+1, aux);
            }
        }
        N_++;
    }
    fclose(f);
    return;
}

void kmap2(Eigen::VectorXd x, float gamma, float lim, char * filename)
{
    int N, i, j;
    int N_ = 1;
    int nz = 0;
    float v, aux, lim_i, lgamma;
    setlocale(LC_NUMERIC, "French_Canada.1252");
    N = x.rows();
    MM_typecode matcode;
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    mm_set_symmetric(&matcode);
    FILE * f = fopen(filename, "w");
    mm_write_banner(f, matcode);

    //std::cout << x << '\n' << std::endl;

    lgamma = log(gamma);
    v = -pow(x(0)-x(1),2) / lgamma;
    lim_i = v * log(lim);
    for(j=0; j<N_; j++){
        if ( pow(x(0)-x(j),2) < -lim_i ){
            nz++;
        }
    }
    N_++;
    for(i=1; i<N; i++){
        if (x(i) < x(i-1)) {
            v = -pow(x(i) - x(i+1),2) / lgamma;
            lim_i = v * log(lim);
        }
        for(j=0; j<N_; j++){
            if ( pow(x(i) - x(j),2) < -lim_i ){
                nz++;
            }
        }
        N_++;
    }
    N_ = 1;
    mm_write_mtx_crd_size(f, N, N, nz);
    v = -pow(x(0)-x(1),2) / lgamma;
    lim_i = v * log(lim);
    for(j=0; j<N_; j++){
        if ( pow(x(0)-x(j),2) < -lim_i ){
            aux = exp(-pow(x(0)-x(j),2)/v);
            fprintf(f, "%d %d %10.3g\n", 1, j+1, aux);
        }
    }
    N_++;
    for(i=1; i<N; i++){
        if (x(i) < x(i-1)) {
            v = -pow(x(i) - x(i+1),2) / lgamma;
            lim_i = v * log(lim);
        }
        for(j=0; j<N_; j++){
            if ( pow(x(i)-x(j),2) < -lim_i ){
                aux = exp(-pow(x(i)-x(j),2)/v);
                fprintf(f, "%d %d %10.3g\n", i+1, j+1, aux);
            }
        }
        N_++;
    }
    fclose(f);
    return;
}
/*
SpMat mmread(char const * filename){
	FILE *fd;
    int M, N, nz;
    MM_typecode matcode;
    fd = fopen(filename, "r");
    mm_read_banner(fd, &matcode);

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    int ret_code = mm_read_mtx_crd_size(fd, &M, &N, &nz);

	std::vector<T> trip;
	trip.reserve(nz);
    int I, J;
    float val;
    for (int i=0; i<nz; i++){
    	fscanf(fd, "%d %d %f\n", &I, &J, &val);
    	trip.push_back(T(I-1,J-1, val));
    	if(I!=J) trip.push_back(T(J-1,I-1, val));
    }
    SpMat m(M, N);
    m.setFromTriplets(trip.begin(), trip.end());
    m.makeCompressed();
    return m;
}
*/
