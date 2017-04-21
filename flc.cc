#include <boost/numpy.hpp>
#include <Eigen/Eigen>
#include "eigen_numpy.h"
#include "mshift.cc"
#include "mm.cc"

#include <iostream>

namespace bp = boost::python;

static const int X = Eigen::Dynamic;

BOOST_PYTHON_MODULE(flc) {
  boost::numpy::initialize();
  SetupEigenConverters();
  bp::def("mshift", mshift< Eigen::Matrix<double,X,X>>);
  bp::def("bmap", bmap); // ok
  bp::def("kmap", kmap); // ok
  bp::def("kmaptau", kmaptau); // ok
  bp::def("kmap2", kmap2);
  bp::def("kmap2tau", kmap2tau);
}
