# Ranking, Order Statistics, and Sorting under CKKS

This repository provides a library for performing ranking, order statistics, and sorting operations using the CKKS (Cheon-Kim-Kim-Song) homomorphic encryption scheme.
Our code is built on top of the OpenFHE library.


## Paper and Citation

This repository accompanies our paper [Efficient Ranking, Order Statistics, and Sorting under CKKS](https://arxiv.org/abs/2412.15126), in publication at USENIX Security '25.

If you use this work, please cite:

```bibtex
@inproceedings{mazzone2025efficient,
  title={Efficient Ranking, Order Statistics, and Sorting under CKKS},
  author={Mazzone, Federico and Everts, Maarten and Hahn, Florian and Peter, Andreas},
  booktitle={34th USENIX Security Symposium (USENIX Security '25)},
  year={2025},
  address={Seattle, WA},
  publisher={USENIX Association},
  month={aug}
}
```


## Contact

For bug reports or inquiries: [contact info](https://people.utwente.nl/f.mazzone).


## Installation (Linux Ubuntu)


### Prerequisites

Install compiler and cmake if needed.

   ```bash
   sudo apt-get install build-essential
   sudo apt-get install cmake
   ```


### OpenFHE Library

To ensure compatibility, we recommend using version **1.1.2** of the OpenFHE library. This version is included within this repository. Follow the steps below to install it:

1. **Navigate to the OpenFHE Directory**

   ```bash
   cd openfhe-development-1.1.2
   ```

2. **Create a Build Directory**

   ```bash
   mkdir build
   ```

3. **Navigate to the Build Directory**

   ```bash
   cd build
   ```

4. **Generate Build Files with CMake**

   ```bash
   cmake ..
   ```

5. **Compile the Library**

   ```bash
   make -j
   ```

6. **Install the Library**

   ```bash
   sudo make install
   ```

7. **Return to the Root Directory**

   ```bash
   cd ../..
   ```

If you do not have sudo access in your machine, specify a different installation path by

   ```bash
   cmake -DCMAKE_INSTALL_PREFIX=~/openfhe ..
   make -j
   make install
   cd ../..
   echo 'export LD_LIBRARY_PATH=~/openfhe/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
   echo 'export LIBRARY_PATH=~/openfhe/lib:$LIBRARY_PATH' >> ~/.bashrc
   echo 'export CMAKE_PREFIX_PATH=~/openfhe:$CMAKE_PREFIX_PATH' >> ~/.bashrc
   echo 'export CPATH=~/openfhe/include:$CPATH' >> ~/.bashrc
   source ~/.bashrc
   ```

For additional support or troubleshooting, refer to the official [OpenFHE Documentation](https://openfhe-development.readthedocs.io/en/latest/).

### Our Library

Follow these steps to compile our library and build the demo and benchmarking executables:

1. **Create a Build Directory**

   Create a directory named `build` in the root of the repository:

   ```bash
   mkdir build
   ```

2. **Navigate to the Build Directory**

   Change into the `build` directory:

   ```bash
   cd build
   ```

3. **Generate Build Files with CMake**

   Use CMake to generate build files based on the configuration specified in `CMakeLists.txt`:

   ```bash
   cmake ..
   ```

4. **Compile the Executables**

   Build the source files to create the executables:

   ```bash
   make -j
   ```

5. **Return to the Root Directory**

   ```bash
   cd ..
   ```

6. **Run the Executables**

   After successful compilation, the following executables will be available in the `build` directory:
   - `demo`
   - `ranking`
   - `minimum`
   - `median`
   - `sorting`

   To run an executable, use the following format:

   ```bash
   ./build/<executable_name>
   ```

   For example, you can run the demo using:

   ```bash
   ./build/demo
   ```

   The `demo` executable demonstrates the functionalities of ranking, finding the minimum, computing the median, and sorting for a given vector.

---



## Benchmarking and Functionality Tests

The executables (`ranking`, `minimum`, `median`, and `sorting`) can be used to benchmark the runtime of the library under various configurations. These executables generate random input vectors, set up the CKKS encryption scheme, and perform the respective operations. By default, the CKKS parameters and approximation degrees are selected automatically, but you may modify them to optimize performance for specific scenarios.

### Usage

Each executable accepts specific arguments as follows:

#### RANKING

```bash
./build/ranking <vector_length> [<tie_correction>] [<single_thread>]
```
- **vector_length**: Length of the input vector (a power of two).
- **tie_correction**: (Optional) Set to `1` to enable tie correction or `0` to disable it (default: `0`).
- **single_thread**: (Optional) Set to `1` to run in single-threaded mode or `0` to use multi-threading (default: `0`).

Example:

```bash
./build/ranking 32 0 1
```

#### MINIMUM

```bash
./build/minimum <vector_length> [<single_thread>]
```
- **vector_length**: Length of the input vector (a power of two).
- **single_thread**: (Optional) Set to `1` to run in single-threaded mode or `0` to use multi-threading (default: `0`).

Example:

```bash
./build/minimum 32 1
```

#### MEDIAN

```bash
./build/median <vector_length> [<single_thread>]
```
- **vector_length**: Length of the input vector (a power of two).
- **single_thread**: (Optional) Set to `1` to run in single-threaded mode or `0` to use multi-threading (default: `0`).

Example:

```bash
./build/median 32 0
```

#### SORTING

```bash
./build/sorting <vector_length> [<single_thread>]
```
- **vector_length**: Length of the input vector (a power of two).
- **single_thread**: (Optional) Set to `1` to run in single-threaded mode or `0` to use multi-threading (default: `0`).

Example:

```bash
./build/sorting 32 1
```

### Automatic Benchmarking

We also provide a shell script `benchmark.sh`, which automatically benchmarks a given functionality for vector length going from 8 to 16384, in both single-threaded and multi-threaded setting. The runtime and memory consumption are stored in a csv file `benchmark.out`, while the log files are stored in the folder `logs`.
```bash
sh ./benchmark.sh <algorithm>
```
- **algorithm**: One of the following: ranking, ranking-tie, minimum, median, sorting.


## Suggested Parameters

The parameter settings used for approximating the comparison and indicator functions in our experiments can be found in the `test-*.cpp` files. Note that:
- The approximation degree of the **comparison function** depends on the precision of the input elements (e.g., sorting 8-bit vs. 16-bit elements).  
- The approximation degree of the **indicator function** depends on the number of input elements (e.g., sorting 32 vs. 64 elements).

We suggest using the **f, g approximation method** with the composition degrees as follows.

**COMPARISON FUNCTION**

| Bit Precision | df | dg |
|--------------|----------|----------|
| 1  | 2 | 1 |
| 2  | 2 | 1 |
| 3  | 2 | 2 |
| 4  | 2 | 2 |
| 5  | 2 | 3 |
| 6  | 2 | 3 |
| 7  | 2 | 4 |
| 8  | 2 | 4 |
| 9  | 2 | 5 |
| 10 | 2 | 5 |
| 11 | 2 | 5 |
| 12 | 2 | 6 |
| 13 | 2 | 6 |
| 14 | 2 | 7 |
| 15 | 2 | 7 |
| 16 | 2 | 8 |

**INDICATOR FUNCTION**
- df = 2
- dg = (log2(vectorLength) + 1) / 2

These settings provide strong correctness guarantees, but they may be overly conservative for many applications.
We encourage you to experiment with lower parameters and/or explore Chebyshev approximation for improved efficiency.  


## Known Issues

**Missing Shared Library**

If you encounter the following error
```error while loading shared libraries: libOPENFHEbinfhe.so.1: cannot open shared object file: No such file or directory```, it means the OpenFHE shared libraries are not in the dynamic linker’s search path.

SOLUTION
1. Temporarily set `LD_LIBRARY_PATH` (valid for the current session only):
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```
2. Make it permanent (applies to future sessions):
```bash
echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```
If OpenFHE is installed in a different location, replace /usr/local/lib with the correct path.




## Docker

We also provide a Docker container, hosted on Docker Hub.
To run the container, follow the instructions at https://hub.docker.com/r/mazzonef/openfhe-statistics .

However, for optimal performance, we recommend using the standard installation.
Running our code in Docker can sometimes introduce runtime overhead, which is particularly noticeable with small vectors.




---
