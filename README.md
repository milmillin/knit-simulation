# Knit Simulation Project

This projects simulates yarn interactions.

## Dependencies and First-Time Setup

The dependencies are STL, GLM (for [upstream](https://github.com/textiles-lab/smobj/) code) , Eigen (for matrix calculations), spdlog (for logging) and [libigl](http://libigl.github.io/libigl/) (for visualization). An optional dependency is [Easy Profiler](https://github.com/yse/easy_profiler).

Since many dependencies are included as git submodules, you need to run `git submodule update --init --recursive` after you clone the repo.

```Bash
cd knit-simulation/
git submodule update --init --recursive
git clone https://github.com/g-truc/glm.git # if you don't have GLM installed before
```

### Easy Profiler
Easy Profiler is used to time the program and help with algorithm optimization. It's intended for developers and not useful for the user.

If you want to avoid this dependency, change `set(ENABLE_EASY_PROFILER TRUE)` to `FALSE` in `CMakeList.txt`.

#### Linux
To insall Easy Profiler, follow the documentation on the [project page](https://github.com/yse/easy_profiler#if-using-cmake), or follow [the following instruction](https://kezunlin.me/post/91b7cf13/) works on Linux:

```Bash
git clone https://github.com/yse/easy_profiler.git

cd easy_profiler && mkdir build && cd build

cmake-gui ..
make -j8
sudo make install
```

#### Windows
Download the release version of easy_profiler [here](https://github.com/yse/easy_profiler/releases) and extract it somewhere. Change the path in `set(EASY_PROFILER_DIR "${PROJECT_SOURCE_DIR}/../easy-profiler")` in `CMakeLists.txt` to match your easy_profiler folder location.

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

    ./knit-simulator <yarns-filename> [-r,--restore] [--no-gui] [-o,--output-dir <output-directory>] [-v,--verbose] [--quiet]

### Options

- `-r,--restore` : (optional) The simulator will not restore the state and history located in the output directory.
- `--no-gui` : (optional) Start the simulator without GUI. (useful when running on a remote machine)
- `-o,--output-dir <output-directory>` : (optional, default = `output/`) Specify the output directory.
- `-v,--verbose` : (optional) Show all logs.
- `--quiet` : (optional) Don't show any logs.

### Output File

The output directory consists of `viewer-state.txt`, `position-xxxxx.yarns`, and `velocity-xxxxx.yarns`.
The `viewer-state.txt` contains number of frames and all the parameters used by the simulator.
This file will be updated at every frame.
