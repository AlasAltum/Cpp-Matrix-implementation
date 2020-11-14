# C-Matrix-implementation
This Class implements C++ matrices to work with 2D matrices, using just 1D arrays. This matrices can be used in CUDA and OpenCL since these frameworks work only with 1D arrays.

## How to use:

Since this class is implemented using templates, the type is going to be settled when you declare a matrix object:

```C++
auto *A = new Matrix<int>(3, 3);
auto *v = new Matrix<int>(3, 1); //vector
auto *u = new Matrix<int>(3); //vector
```

In order to use this, you just need to include this file:

```C++
#include "matrix.h"
```

In main.cpp you can find some testing of this library.

# Remarks

This was done as an introductory class to C++ for a course.
