# Master Thesis

##Incompact3d Compiling
To compile and run the incompact3d:

```bash
git clone https://github.com/wwdhd/master_thesis.git
cd Incompact3d_osc180
cd ../
export FC=mpif90
cmake -S . -B build
cd build
cmake --build . -j {number of cores}
cd ../examples/Channel/

