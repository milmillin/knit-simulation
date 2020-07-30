# Knit Simulation Project

This projects simulates yarn interactions.

## Dependencies and First-Time Setup

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`. An optional dependency is [Easy Profiler](https://github.com/yse/easy_profiler).

For now, we also require [glm](https://github.com/g-truc/glm).

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

```Bash
cd knit-simulation/
git clone https://github.com/libigl/libigl.git
git clone https://github.com/g-truc/glm.git
```

Easy Profiler is used to time the program and help with algorithm optimization. It's intended for developers and not useful for the user. If you want to avoid this dependency, comment out `#define ENABLE_EASY_PROFILER` in `main.cpp`.

To insall Easy Profiler, follow the documentation on the [project page](https://github.com/yse/easy_profiler#if-using-cmake).

[This instruction](https://kezunlin.me/post/91b7cf13/) works on Linux:

```Bash
git clone https://github.com/yse/easy_profiler.git

cd easy_profiler && mkdir build && cd build

cmake-gui ..
make -j8
sudo make install
```

## Compile

Compile this project using the standard cmake routine:

```Bash
mkdir build
cd build
cmake ..
cmake --build .
```

This should find and build the dependencies and create a `knit-simulator` binary.

## Run

From within the `build` directory just issue:

    ./knit-simulator

A glfw app should launch displaying a 3D cube.

## Note

You might want to change the path to example yarns file in `main.cpp` according to your build directory.

    yarn = file_format::Yarns::load("../../../helloworld.yarns");
