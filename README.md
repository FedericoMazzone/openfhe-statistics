# Ranking, Order Statistics, and Sorting under CKKS

This repository provides a library for performing ranking, order statistics, and sorting operations using the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption scheme.

## Installation

To install OpenFHE, the required homomorphic encryption library, please refer to the documentation available at [OpenFHE Documentation](https://openfhe-development.readthedocs.io/en/latest/).

## Compilation Instructions

Follow these steps to compile the library and run the provided demo:

1. **Set Files to Compile:** Modify the `CMakeLists.txt` file to specify the source files you want to compile. By default, a `demo.cpp` file is included as an example.

2. **Create Build Directory:** Create a directory named `build` in the root of the repository.

    ```bash
    mkdir build
    ```

3. **Navigate to Build Directory:** Move into the `build` directory.

    ```bash
    cd build
    ```

4. **Run CMake:** Run CMake to generate build files based on the configurations specified in `CMakeLists.txt`.

    ```bash
    cmake ..
    ```

5. **Build Executable:** Compile the source files and build the executable.

    ```bash
    make
    ```

6. **Run Program:** Once the compilation is successful, the executable named `program` will be available in the `build` directory. You can run it to execute the provided demo.
