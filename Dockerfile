FROM ubuntu:22.04

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libgmp-dev \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy files
COPY data /app/data
COPY openfhe-development-1.1.2 /app/openfhe-development-1.1.2
COPY src /app/src
COPY benchmark.sh /app/benchmark.sh
COPY CMakeLists.txt /app/CMakeLists.txt
COPY LICENSE /app/LICENSE
COPY PreLoad.cmake /app/PreLoad.cmake
COPY README.md /app/README.md

# Build OpenFHE
WORKDIR /app/openfhe-development-1.1.2
RUN mkdir build && cd build && \
    cmake .. && \
    make -j && \
    make install

# Set library path
ENV LD_LIBRARY_PATH=/usr/local/lib

# Build the project
WORKDIR /app
RUN mkdir build && cd build && \
    cmake .. && \
    make -j

# Set default command
CMD ["/bin/bash"]
