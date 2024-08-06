# Master Thesis

To run the incompact3d:

`git clone https://github.com/wwdhd/master_thesis.git`

::cd Incompact3d_osc180::

cd src #if needed

cd ../

export FC=mpif90

cmake -S . -B build

cd build

cmake --build . -j 8

cd ../examples/Channel/
