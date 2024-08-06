# Master Thesis

## Incompact3d Compiling

To compile and run the incompact3d (substitute ncores with the number of desired cores for compiling):

```bash
git clone https://github.com/wwdhd/master_thesis.git
cd Incompact3d_osc180
cd ../
export FC=mpif90
cmake -S . -B build
cd build
cmake --build . -j {ncores}
cd ../examples/Channel/
```
If you are using HPC environment, loading cmake module is important. To run the channel flow code in the local computer:

```bash
mpirun -np {ncores} ../../build/bin/xcompact3d > x3d2.log &
tail -f x3d2.log
```

