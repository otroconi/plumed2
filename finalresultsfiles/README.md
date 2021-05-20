The files present in the current folder corresponds to pure SimpleMD.cpp implementation, i.e. in the file plumed2/src/cltools/plumed.dat comment with # all lines INCLUDE FILE=... and write the corresponding name of the COLVAR file. For example, if after the pure Lennard Jones potential you add the bonds effects in the SimpleMD.cpp, then you modify the COLVAR file name in the plumed.dat file by COLVAR-B and run:

```
  make obj
  cd ../lib
  make lib plumed
  cd -
  plumed simplemd < input
```

Hence, in the COLVAR-B file you will see the effects of the Lennard Jones and the bonds interactions.

Additionally, you add the angles effects, follow the steps previously described, and you will see the file COLVAR-A with the effects of Lennard Jones, bonds and angles interactions. 

Once you add an interaction you can compare the COLVAR files with those present in the folder plumed/runfolder, for example after you add the bonds interaction you can compare your COLVAR-B file with the file runfolder/runbonds/COLVAR-B that corresponds to the plumed result, i.e. that file corresponds to run the original SimpleMD.cpp code with ONLY the Lennard Jones interaction and including in the plumed.dat file the INCLUDE FILE=bonds.dat. 
