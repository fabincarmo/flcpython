#include <iomanip>
#include "mmio.h"
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <tuple>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;

void bmap(Eigen::VectorXd v, char * filename);
void kmap(Eigen::VectorXd x, float gamma, float lim, char * filename);
void kmap2(Eigen::VectorXd x, float gamma, float lim, char * filename);
//SpMat mmread(char const * filename);
