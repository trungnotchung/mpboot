# Compilation guide
## MPBoot
### Downloading source code
You can clone the source code from GitHub with:

`git clone https://github.com/diepthihoang/mpboot.git`
### Compiling under Linux
1. Create folder **build** outside folder **mpboot**.
2. Open Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.  
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy**, run command:  
`./mpboot-avx -s example.phy`

### Compiling under Mac OS X
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.  
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy**, run command:  
`./mpboot-avx -s example.phy`

### Compiling under Windows
* Requirements:  
  * cmake version >= 3.21
  * TDM-GCC
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake -G "MinGW Makefiles" -DIQTREE_FLAGS=avx ../mpboot`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture.  
> Due to having conflicts with **Vectorization**, please not using **Clang** to configure source code.

5. Run command `mingw32-make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate two executable files named **mpboot-avx** and **mpboot-avx-click**
> In case of running MPBoot on SSE architecture, the executable files are named **mpboot** and **mpboot-click**.

6. To analyst file **example.phy**:
* Run command `./mpboot-avx -s example.phy`
* Or:
    1. Open **mpboot-click** by double click on it
    2. Press **y** to start enter command
    3. Type **-s** and press **Enter**
    4. Press **e** to continue entering command
    5. Type **example.phy** and press **Enter**
    6. Press **y** to finish command


  
## MPBoot-MPI
### Downloading source code
You can clone the source code from GitHub with:

`git clone https://github.com/diepthihoang/mpboot.git`  
Switch to branch which contains MPBoot-MPI source code:
* `git checkout mpboot-mpi-sync` for synchronous version
* `git checkout mpboot-mpi-async` for asynchronous version

### Compiling under Linux
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../source -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx**
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy** with 2 processes, run command:  
`mpirun -np 2 ./mpboot-avx -s example.phy`
> Option **np** specifies the number of processes used to run MPBoot-MPI

### Compiling under Mac OS X
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake ../mpboot -DIQTREE_FLAGS=avx -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx`
> Replace avx by sse4 in above command if you decide to run MPBoot on SSE architecture

5. Run command `make -j4` to compile source code with 4 processes. Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx**
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.

6. To analyst file **example.phy** with 2 processes, run command:  
`mpirun -np 2 ./mpboot-avx -s example.phy` 
> Option **np** specifies the number of processes used to run MPBoot-MPI

### Compiling under Windows
* Requirements:  
  * cmake version >= 3.21
  * TDM-GCC
  * MSMPI
1. Create folder **build** outside folder **mpboot**.
2. Open a Terminal.
3. Change directory to **build**
4. Configure source code with CMake:  
`cmake -G "MinGW Makefiles" -DIQTREE_FLAGS=mpiavx ../mpboot`
> Replace mpiavx by mpisse4 in above command if you decide to run MPBoot on SSE architecture.  
> Due to having conflicts with **Vectorization**, please not using **Clang** to configure source code.

5. Run command `mingw32-make -j4` to compile source code with 4 processes:  
> Option **j** specifies the number of processes used to compile source code with **make**.
  
The compiler will generate an executable file named **mpboot-avx** 
> In case of running MPBoot on SSE architecture, the executable file is named **mpboot**.
  

6. To analyst file **example.phy** with 2 processes, run command:  
`mpiexec -n 2 ./mpboot-avx -s example.phy`
> Option **n** specifies the number of processes used to run MPBoot-MPI

## MPBoot-Placement

MPBoot-Placement is a specialized tool for phylogenetic placement that efficiently places new samples onto an existing phylogenetic tree using parsimony-based algorithms. It processes VCF files and can handle large datasets through intelligent batch processing.

### Downloading source code
You can clone the source code from GitHub with:

`git clone https://github.com/diepthihoang/mpboot.git`

### Compilation
MPBoot-Placement uses the **same compilation process** as standard MPBoot (see sections above). The executable generated is the same `mpboot-avx` (or `mpboot` for SSE).

### Usage

**Basic placement command:**
```bash
./mpboot-avx -pp_on -s variants.vcf -pp_tree reference.newick -pp_n 1000 -pp_k 10
```

**Placement-specific options:**
- `-pp_on` - Enable phylogenetic placement mode (REQUIRED)
- `-pp_tree <treefile>` - Reference phylogenetic tree in Newick format (REQUIRED)
- `-pp_n <number>` - Number of existing sequences already in the tree
- `-pp_k <number>` - Number of new sequences to place onto the tree
- `-pp_spr` - Enable SPR (Subtree Pruning and Regrafting) optimizations
- `-pp_test_optimize` - Verify that the original tree structure is preserved

**Input requirements:**
- `-s <vcf_file>` - Input VCF file containing genetic variants (REQUIRED)
- VCF file must have standard format with proper header lines
- Reference tree in Newick format containing existing sequences

**Example commands:**

1. **Basic placement:**
   ```bash
   ./mpboot-avx -pp_on -s variants.vcf -pp_tree reference_tree.newick -pp_n 500 -pp_k 50
   ```

2. **With SPR optimization:**
   ```bash
   ./mpboot-avx -pp_on -s variants.vcf -pp_tree reference_tree.newick -pp_n 500 -pp_k 50 -pp_spr
   ```

3. **With verification (testing mode):**
   ```bash
   ./mpboot-avx -pp_on -s variants.vcf -pp_tree reference_tree.newick -pp_n 500 -pp_k 50 -pp_spr -pp_test_optimize
   ```