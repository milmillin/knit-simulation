# Knit Simulation Project

This projects simulates yarn interactions.

## Dependencies and First-Time Setup

The dependencies are STL, GLM (for [upstream](https://github.com/textiles-lab/smobj/) code) , Eigen (for matrix calculations), spdlog (for logging) and [libigl](http://libigl.github.io/libigl/) (for visualization). An optional dependency is [Easy Profiler](https://github.com/yse/easy_profiler).

The cmake build system will attempt to find libigl according to environment variables (e.g., `LIBIGL`) and searching in common desitinations (e.g., `/usr/local/libigl/`). If you haven't installed libigl before, we recommend you to clone a copy of libigl right here:

```Bash
cd knit-simulation/
git clone https://github.com/libigl/libigl.git
git clone https://github.com/g-truc/glm.git
```

### Easy Profiler

Easy Profiler is used to time the program and help with algorithm optimization. It's intended for developers and not useful for the user.

If you want to avoid this dependency, comment out `#define USE_EASY_PROFILER` in `easy_profiler_stub.h`. In `CMakeList.txt`, remove `find_package(easy_profiler REQUIRED)` and remove `easy_profiler` in `target_link_libraries`.

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

## Usage

    ./knit-simulator yarns-filename [--no-restore] [--no-gui] [--output output-directory]

### Options

- `--no-restore` : (optional) The simulator will not restore the state and history located in the output directory.
- `--no-gui` : (optional) Start the simulator without GUI. (useful when running on a remote machine)
- `--output output-directory` : (optional, default = `output/`) Specify the output directory.

### Output File

The output directory consists of `viewer-state.txt`, `position-xxxxx.yarns`, and `velocity-xxxxx.yarns`.
The `viewer-state.txt` contains number of frames and all the parameters used by the simulator.
This file will be updated at every frame.
