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


  
## Placement Optimizer
MPBoot includes a parsimony-based **placement + post-placement optimization** pipeline.

Given an existing reference tree and a set of new samples (VCF), it:
1. **Places** the new samples onto the reference tree (online placement using a mutation-annotated tree).
2. **Optimizes** the placed tree with three phases:
   - **SPR** — exponential-radius hill-climbing (`r = 1, 2, 4, …, max_radius`) with an O(M̄) sparse delta evaluator, batched move application, and verify-by-recompute rollback.
   - **Ratchet** *(optional)* — Nixon parsimony ratchet (pattern reweighting ×2 on ~25% of patterns) combined with Metropolis acceptance (`exp(-Δ·β)`, β=2) for escaping local minima.
   - **TBR** *(optional)* — interleaved TBR rounds (1 TBR + 10 SPR passes) for moves SPR cannot reach.

All three phases share a byte-per-pattern Fitch engine with per-vertex diff vectors and an O(1) LCA table.

### Build
Same as MPBoot above — no separate target. The placement-optimizer code is compiled into the same `mpboot-avx` (or `mpboot` on SSE) binary.

### Inputs
| Input | Flag | Format |
|---|---|---|
| Reference tree | `-pp_tree <file>` | Newick (`.treefile`) |
| Reference alignment + new samples | `-s <file>` | VCF |

### Run — placement only
```
./mpboot-avx -s <vcf_file> -pp_tree <reference_tree> -pp_on -pp_k <num_new_samples> -pp_n <num_existing_samples>
```
Runs placement, does NOT run optimization.

### Run — placement + full optimization (SPR + Ratchet + TBR)
```
./mpboot-avx -s <vcf_file> -pp_tree <reference_tree> \
  -pp_on -pp_k <num_new_samples> -pp_n <num_existing_samples> \
  -pp_optimize \
  -pp_max_radius 32 \
  -pp_ratchet_iters 30 -pp_ratchet_runs 3 -pp_ratchet_seed 42 \
  -pp_tbr_iter 5 -pp_tbr_max_radius 5 \
  -pp_wall_seconds 300
```
Each phase activates only if its iteration count is nonzero — set `-pp_ratchet_iters 0` or `-pp_tbr_iter 0` to skip a phase.

### Key flags
| Flag | Default | Purpose |
|---|---|---|
| `-pp_on` | off | Enable the placement pipeline |
| `-pp_tree <file>` | — | Reference Newick tree |
| `-pp_k <n>` | — | Number of new samples |
| `-pp_n <n>` | — | Total of existing samples |
| `-pp_optimize` | off | Enable post-placement optimization |
| `-pp_max_radius <r>` | 32 | Max SPR radius (`0` = unbounded) |
| `-pp_wall_seconds <s>` | 0 | Wall-clock cap for the whole optimization (`<=0` = no cap) — checked at pass / round / DFS-inner levels |
| `-pp_ratchet_iters <n>` | 0 | Ratchet iterations per run (`0` = disabled) |
| `-pp_ratchet_runs <K>` | 1 | Independent ratchet restarts; the best-of-K tree is kept |
| `-pp_ratchet_seed <s>` | 42 | RNG seed for pattern reweighting |
| `-pp_tbr_iter <n>` | 0 | Outer TBR rounds (`0` = disabled) |
| `-pp_tbr_max_radius <r>` | 5 | BFS depth from the bisection scar when searching reconnection edges |

  
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
