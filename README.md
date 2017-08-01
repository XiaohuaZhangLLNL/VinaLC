# MPI and Multi-threading mix Parallel AutoDock Vina.


## 1. Compile the program.

This source code is configured to run on LLNL LC machines. 

### 1.1 Boost library (www.boost.org) is require for Vina.
to install Boost library please follow the step in the Boost document. 

```
./bootstrap.sh --prefix=path/to/installation/prefix
./b2 install
```

Beside the standard installation, Boost MPI binding also need to be turn on.
copy tools/build/v2/user-config.jam to your home directory. In the file 
specify the mpi compiler you want to use:

```
using mpi : /usr/local/tools/mvapich-gnu/bin/mpicxx ;
```

On quartz:
```
module load boost/1.62.0
```
### 1.2 Obtain the code 

The code can be download from:

https://lc.llnl.gov/bitbucket/projects/XZR/repos/vinalc/browse

by git:
```
git clone ssh://git@cz-bitbucket.llnl.gov:7999/xzr/vinalc.git
```


### 1.3 Installation

Installing it by using cmake is straight forward:

```
cd vinalc
cmake . -DCMAKE_INSTALL_PREFIX:PATH=/usr/gapps/bbs/vinalc-quartz
make
make install
```

## 2. Input files.
In the vina-release/build/linux/debug, there is a small test case contains input files:

### 2.1 receptor

file contains list of receptor file name:  recList.txt
two receptor pdbqt files: 

```
1KIJ_protH.pdbqt 1KIJ_protH1.pdbqt
```

### 2.2 ligand

file contains list of ligand file name: ligList.txt
a "data" directory contains two ligand pdbqt files:

```
data/
ligands1.pdbqt  ligands2.pdbqt
```

### 2.3 geometry file to save the docking grid information.
each line is coresponding to each receptor in recList.txt.
each line has six number where first 3 are center of active site and last 3 are grids dimension.
the file is arranged as:
```
x_center  y_center  z_center  x_grid  y_grid  z_grid
x_center  y_center  z_center  x_grid  y_grid  z_grid
....
```

### 2.4 The pdbqt files for receptors and ligands are prepared from their pdb file by mgltools
(http://mgltools.scripps.edu/)
Two python scripts in MGLTools-1.5.6rc2/MGLToolsPckgs/AutoDockTools/Utilities24/ are used

for receptor
```
prepare_receptor4.py -r rec.pdb -o rec.pdbqt -A checkhydrogens
```

for ligand
```
prepare_ligand4.py -l ligand.pdb
```

## 3. Running program

### 3.1 To run the program with slurm in debug mode:

```
srun -N4 -n4 -c12 -ppdebug ./vina --recList recList.txt --ligList ligList.txt --geoList geoList.txt

-N4: 4 nodes will use
-n4: 4 tasks will  each task run on one node
-c12: 12 threads will run on each node 
-ppdebug: use debug mode
```

### 3.2 Vina program option:
```
[zhang30@quartz1538 bin]$ ./vinalc --help

Input:
  --recList arg               receptor list file
  --fleList arg               flex part receptor list file
  --ligList arg               ligand list file
  --geoList arg               receptor geometry file
  --exhaustiveness arg (=8)   exhaustiveness (default value 8) of the global 
                              search (roughly proportional to time): 1+
  --granularity arg (=0.375)  the granularity of grids (default value 0.375)
  --num_modes arg (=9)        maximum number (default value 9) of binding modes
                              to generate
  --seed arg                  explicit random seed
  --randomize                 Use different random seeds for complex
  --energy_range arg (=2)     maximum energy difference (default value 2.0) 
                              between the best binding mode and the worst one 
                              displayed (kcal/mol)
  --useScoreCF                Use score cutoff to save ligand with top score 
                              higher than certain critical value
  --scoreCF arg (=-8)         Score cutoff to save ligand with top score higher
                              than certain value (default -8.0)

Information (optional):
  --help                      display usage summary

Error: Total process less than 2

```

### 3.3 Run with different options

```
srun -N4 -n4 -c12 ./vina --recList recList.txt --ligList ligList.txt --geoList geoList.txt --exhaustiveness 12
srun -N4 -n4 -c12 ./vina --recList recList.txt --ligList ligList.txt --geoList geoList.txt --exhaustiveness 12 --granularity 0.333
...

options, --recList --ligList  --geoList must be specified  
```

## 4. The programs to prepare the input PDBQT file.

### 4.1 babel
If you start with SDF files for the ligands, you can use the babel to convert the SDF to PDBQT directly:
```
babel -isdf <ligand-file-name>.sdf -opdbqt <ligand-file-name>.pdbqt
```

### 4.2 MGLTools(http://mgltools.scripps.edu/)
The AutoDock developer team provides graphic user interface, AutoDockTools (ADT), to prepare the input files. The receptor input file MUST use ADT to convert the file format.

There is a Vina video tutorial to show how to use ADT to prepare receptor, ligand, and determine the grid size that use in the program.

http://vina.scripps.edu/tutorial.html

An important thing to remember when calculate the grid size: 
```
x_grid=<number of poiont in x-dimension>*spacing
y_grid=<number of poiont in y-dimension>*spacing
z_grid=<number of poiont in z-dimension>*spacing
```

spacing in ADT is equal to granularity in Vina.



