# Knit Simulation Project

This is a Knit Simulation Project.

## Dependencies and First-Time Setup

The only dependencies are stl, eigen, [libigl](http://libigl.github.io/libigl/) and
the dependencies of the `igl::opengl::glfw::Viewer`.

For now, we also require [glm](https://github.com/g-truc/glm).

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

    cd knit-simulation/
    git clone https://github.com/libigl/libigl.git
    git clone https://github.com/g-truc/glm.git

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    cmake --build .

This should find and build the dependencies and create a `example_bin` binary.

## Run

From within the `build` directory just issue:

    ./knit-simulator

A glfw app should launch displaying a 3D cube.
