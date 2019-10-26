# N-Body

 Authors:
     - Weizhi Xu   (weizhix)  752454
     - Zijun Chen  (zijunc3)  813190

## File Structure
```
/
├── src/             // main codes
├── slurm/           // slrum file for running on Spartan
├── results/         // raw results from experiments
├── tools/
│   ├── gui.py       // gui for playing simulation results
│   └── evaluate.py  // tool for computing the accuracy of results
└── generator.py     // test data generator
```


## Usage

Build none CUDA version
```
make nbody
```

Build CUDA version
```
make nbody-with-cuda
```

You need to use make clean when switching targets.

```
Usage: nbody <algorithm name> <num of steps> <time of each step> <full output> <input_file> <output file>
Example: nbody seq_naive 100 0.01 1 input.in output.out
Supported algorithms:
	seq_naive
	seq_quad_tree
	mpi_openmp_naive
	mpi_quad_tree
	cuda_single_naive
	cuda_mpi_naive
```

