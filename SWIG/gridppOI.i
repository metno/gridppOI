%module gridppOI
%include "typemaps.i"
%include "std_vector.i"
namespace std {
  %template(IntVector) vector<int>;
  %template(DoubleVector) vector<double>;
  %template(FloatVector) vector<float>;
  %template(FloatVector2) vector<vector<float> >;
}
%apply std::vector<std::vector<float> >& OUTPUT { std::vector<std::vector<float> >& output };
%{
/*  Put header files here or function declarations like below */
#include "gridppOI.h"
%}
%include "gridppOI.h"
