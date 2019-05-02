// STL
#include <cassert>

// Pybind11
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

// Local
#include "BufferWrapper.h"
//#include "Class2Wrapper.h"
//#include "Class3Wrapper.h"

namespace py = pybind11;

PYBIND11_MODULE(pysiril, m) {
  m.doc() = "pysiril : pybind11 siril binding";

  py::class_<BufferWrapper<float>> Buffer(m, "Buffer",
    py::dynamic_attr());
  Buffer
    .def(py::init<int>())
    .def("initialize", &BufferWrapper<float>::initialize)
    .def("set_image", &BufferWrapper<float>::set_image);

