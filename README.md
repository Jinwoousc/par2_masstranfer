# PAR² MassTransfer
PAR² is a Lagrangian solute transport simulator using a parallelized Random Walk Particle Tracking (RWPT) method. This has been updated to incorporate masstransfer.

# What is in the MassTransfer update?
Mass Transfer Models (e.g., partitioning or one-site models) is incorporated for each grid block to simulate more realistic solute transport when mass transfer between mobile and immobile phases occurs.

For instance, the figures below clearly show the results of mass transfer. The first (upper left) and second (upper right) figures present solute transport without mass transfer in homogeneous and heterogeneous porous media, respectively. The third (lower left) and fourth (lower right) figures exhibit solute transport with mass transfer (one-site model) in the same homogeneous and heterogeneous porous media. You can see the immobile phase (brown dots) left along the way of the mobile phase (blue dots).

![no_masstransfer_homogeneous_case](/README_images/no_masstransfer_homogeneous_case.png) ![no_masstransfer_heterogeneous_case](/README_images/no_masstransfer_heterogeneous_case.png)
![one-site_model_homogeneous_case](/README_images/one-site_model_homogeneous_case.png) ![one-site_model_heterogeneous_case](/README_images/one-site_model_heterogeneous_case.png)

## Getting Started
You can download the latest release of PAR² for Windows executable from [here](https://github.com/Jinwoousc/par2_updated/tree/main/Build/Release). Make sure your computer is equipped with an NVIDIA GPU and NVIDIA drivers are updated. You can control the simulation parameters through a YAML configuration file. Look inside the Examples folder to get started.

## Build
The following software and libraries must be installed:

* cmake (version 3.8 or higher)
* CUDA Toolkit (tested with version 9.0)
* yaml-cpp library
* spdlog library (included in the source code)

### Compile on Linux
1. Make sure to have a valid c++ compiler (e.g., gcc)
2. Create a build directory:

        mkdir Build  
        cd Build

2. Create makefile:

        cmake -DCMAKE_BUILD_TYPE=Release ..

3. Compile:

        make

### Compile on Windows
1. Make sure to have Visual Studio 2017 installed
2. Create a build directory:

        mkdir Build  
        cd Build

2. Create MSVC solution:

        cmake -G "Visual Studio 15 2017 Win64" -T v140 -DCMAKE_BUILD_TYPE=Release -DYAML_ROOT=C:/path/to/yaml-cpp ..

3. Compile using the Developer Command Prompt for VS:

        devenv /Build Release par2.sln

## Citations
Rizzo, C. B., Nakano, A., & de Barros, F. P. J. (2019). [PAR²: Parallel Random Walk Particle Tracking Method for solute transport in porous media](https://doi.org/10.1016/j.cpc.2019.01.013) Computer Physics Communications, 239, 265-271.

Salamon, P., Fernàndez‐Garcia, D., & Gómez‐Hernández, J. J. (2006). [Modeling mass transfer processes using random walk particle tracking](https://doi.org/10.1029/2006WR004927) Water Resources Research, 42(11).
