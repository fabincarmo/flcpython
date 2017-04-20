#include <iostream>
#include <eigen3/Eigen/Dense>

using namespace Eigen;
using namespace std;

template <typename M> M mshift(M in, int down, string dir )
{
  M out =  MatrixXd::Ones(in.rows(), in.cols()) * 1e-16;
  if (dir=="ud"){
    int d;
    if (!down) return in;
    if (down > 0) d = down % in.rows();
    else d = -down;
    int rest = in.rows() - d;
    if (down > 0) out.topRows(rest) = in.bottomRows(rest);
    else out.bottomRows(rest) = in.topRows(rest);
    return out;
  }
  if (dir=="lr"){
    int d;
    if (!down) return in;
    if (down > 0) d = down % in.cols();
    else d = -down;
    int rest = in.cols() - d;
    if (down > 0) out.rightCols(rest) = in.leftCols(rest);
    else out.leftCols(rest) = in.rightCols(rest);
    return out;
  }
}
