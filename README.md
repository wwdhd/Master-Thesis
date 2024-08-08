# Master Thesis

## Incompact3d Compiling
Before compiling, cmake and mpich modules must be installed into Ubuntu, using `sudo apt install mpich` and `sudo apt install cmake`. If they are unable to be downloaded, make sure Ubuntu is already up to date using `sudo apt update` and then `sudo apt upgrade`. If you are using HPC environment, loading cmake and mpich module (or openMPI in some cases) is important.

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

 To run the channel flow code in the local computer:

```bash
mpirun -np {ncores} ../../build/bin/xcompact3d > x3d2.log &
tail -f x3d2.log
```
In HPC environment, the same command is used, just ensure the job scheduler format is followed.

## Py4Incompact3d Download and Compiling
Before compiling, it is important to have the Python module, downloaded using ` sudo apt install python3-pip`. To compile and install Py4Incompact3d in Ubuntu:
```bash
git clone https://github.com/xcompact3d/Py4Incompact3D
cd Py4Incompact3D
pip install .
pip install scipy
```

