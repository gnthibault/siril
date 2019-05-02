// STL
#include <cassert>

// Local
//#include "Buffer.h"

namespace py = pybind11;

template<typename T>
class BufferWrapper {
 public:
  BufferWrapper(int) {};
  void initialize(py::array_t<T> image) {
    auto buffer = image.request();
    if (buffer.ndim < 2) or (buffer.ndim > 3) {
      throw std::runtime_error("BufferWrapper::initialize : "
        "Number of dimensions must be 2 or 3");
    }
    T* ptr = static_cast<T *>(buffer.ptr);
    m_pData = ptr;
    m_dataSize = buffer.size;
    //buffer.shape[0]*buffer.shape[1]);
  }
  void set_image(py::array_t<T> image) {
    auto buffer = image.request();
    if (buffer.ndim < 2) or (buffer.ndim > 3) {
      throw std::runtime_error("BufferWrapper::set_image : "
        "Number of dimensions must be 2 or 3");
    }
    m_pData = ptr;
    m_dataSize = buffer.size;
  }

 protected:
  T* m_pData;
  size_t m_dataSize;
};
