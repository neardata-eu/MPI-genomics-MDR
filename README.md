# MPI-genomics-MDR

## Project Description

MPI-genomics-MDR is an application for identifying genomic variant interactions related to complex diseases. It utilizes Multifactor Dimensionality Reduction (MDR), a non-parametric statistical method, to detect and characterize nonlinear interactions.  Utilizing Message Passing Interface (MPI) functions, MPI-genomics-MDR efficiently distributes computational tasks across multiple nodes, harnessing the collective computing resources for enhanced scalability and performance.

## Prerequisites

Before compiling and running the program, ensure that you have the following installed:

- MPI implementation (e.g., Open MPI, MPICH)
- C compiler (e.g., GCC)
- Set the `ROOT_PATH` in the `scripts/mdr-mpi.h` file

## Usage

To compile the `mdr-mpi` program, follow these steps:

1. Clone this repository to your local machine:
   ```
   git clone https://github.com/neardata-eu/MPI-genomics-MDR.git
   ```
2. Navigate to the scripts directory:
   ```
   cd repository/scripts
   ```

3. Type `make` to compile:
   ```
   make
   ```

## Running the Experiments

To run the experiments:

1. Execute the compiled program using MPI:
   ```
   mpiexec -n <num_processes> ./mdr-mpi
   ```
   Replace `<num_processes>` with the desired number of MPI processes

## License

Copyright 2024 Andrés Benavides A and Gonzalo Gómez-Sánchez

Licensed under the Apache License, Version 2.0 (the "License"). You may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0).

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
