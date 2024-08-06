# Master Thesis

To run the incompact3d:

|git clone https://github.com/wwdhd/master_thesis.git

| *Yuku haru ya*
| *tori naki uo no*
| *me wa namida*
| -- **Matsuo Bashō**, The Narrow Road to Oku (*Oku no Hosomichi*),
 Tokyo, 1996, p. 23 (Translation: Donald Keene)
| Spring is passing by!
| Birds are weeping and the eyes
| of fish fill with tears.

`cd Incompact3d_osc180`

`cd ../`

`export FC=mpif90`

`cmake -S . -B build`

`cd build`

`cmake --build . -j 8`

`cd ../examples/Channel/`
