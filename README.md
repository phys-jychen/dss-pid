# DSS-PID, a Simple Tool for PID of DarkSHINE Baseline ECAL Design

## Author
Ji-Yuan CHEN (SJTU; <jy_chen@sjtu.edu.cn>)

## Description
This program is designed for the PID of DarkSHINE Baseline ECAL Design. By reconstructing variables describing the topology of the hadronic and electromagnetic showers, we can perform PID with the help of BDT, using the TMVA (the **T**oolkit for **M**ulti-**V**ariate data **A**nalysis with ROOT) package.

For a more intuitive understanding of the shower topology, a python program is also included for event display. You can also use `DDis` in DarkSHINE simulation framework.

## Usage
At first, you can run
```shell
iHit -h[elp]
# or
iRec -h[elp]
# or
iBDT -h[elp]
```
to display help information of analysis, and
```shell
python display.py -h
```
for that of event display (mind the order of `display.py` and `-h`!). For more detail, please refer to the following instructions. :stuck_out_tongue:

### Sample Preparation
After generating MC samples with `DSimu` in the DarkSHINE simulation framework, assign a large number (e.g. 999999) to the option `RecECAL.Advance` in the configuration file for `DAna`. Then run `DAna`, and the original hit information will be stored in the branch `ECAL_ECell_XYZ`.

### Selecting Hits
The desired hit information are in the form of `Hit_X`, `Hit_Y`, etc., which should be converted from `ECAL_ECell_XYZ`. To achieve this, run
```shell
iHit -s -f [file]
```
if the tree in the original ROOT file has default name `dp`. Otherwise, to specify other tree names, run
```shell
iHit -s -f [file] -t [tree]
```

After this, an output file whose name has a prefix “sel” is created in your current directory, regardless of where the original ROOT file is. Notice that this step keeps only the branch `ECAL_ECell_XYZ`, since some other branches lead to failure in collecting original hits!

Next, it is time to store explicit hit information:
```shell
iHit -f [file]
# or
iHit -f [file] -t [tree]
```

Then, an output file whose name has a prefix “sel” is created in your current directory. Other branches have already been deleted to save space.

### Adding Reconstructed Variables
In any directory, execute:
```shell
iRec -f [file]
# or
iRec -f [file] -t [tree]
```

The name of the output file is given a prefix “rec” in your current directory, and the original branches are not kept for the sake of saving space. If you need to add some new variables or modify the definitions of some of them, please go to the file `src/Variables.cxx`.

### Performing BDT
Before you begin, make sure that the variables as well as the ROOT files listed in `bdt.cxx` are all present (you can also modify this file to meet your own needs). Then, execute:
```shell
iBDT -r
# or
iBDT -r -t [tree]
```

Eventually, you can see the output file containing the data obtained during training and test (`TMVAMulticlass.root`) in your current directory.

Possibly you need to know the performance of BDT on the validation dataset. In this case, execute:
```shell
iBDT -v -f [file]
# or
iBDT -v -f [file] -t [tree]
```

Then, the BDT response is stored in the output ROOT file, whose name has a prefix “bdt”, in your current directory. While modifying `src/BDT.cxx`, make sure that the input variables are identical to those in `bdt.cxx`, including the order!

### Event Display
In the directory you have installed, run
```shell
./event_display.sh
```
to obtain a figure of event display, which will be saved in a directory assigned in this shell script (default: `figs/`).

You can modify `event_display.sh` to meet your own needs: input ROOT file, tree name, title of the figure, ID of the event, directory to save the output file, name of the output file, and instantly show the figure or not.

## Environment Set-up
This project requires CMake version >= 3.17. If you are working on the cluster of INPAC, IHEP, etc., the environment can be easily set up by simply executing
```shell
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh
```
(This command has been included in `setup.sh`.)

Then, the environment with CMake 3.20.0 and ROOT 6.26/04 is configured. As long as neither compilation errors are raised, nor the CMake version requirement is met, other versions on the LCG are also acceptable. :stuck_out_tongue:

### Notice
On WSL2, the above environment cannot be used for event display, since error occurs with X11! Instead, the environment used while designing this program was with
- Python: 3.11.5
- Matplotlib: 3.7.2
- NumPy: 1.24.3
- UpROOT: 5.1.1

Other environments have not been fully tested.

## Installation & Compilation
Having set up the environment, this project can be cloned from GitHub and compiled as usual:
```shell
git clone git@github.com:phys-jychen/dss-pid.git
cd dss-pid
mkdir build
cd build
cmake ..
make -j100    # Just do it!
source setup.sh
```

Every time you log in to the cluster, and right before the first time of running this program, remember to execute
```shell
source <build_dir>/setup.sh
```

By now, the compilation have been finished. Prepare your datasets, and have fun! :relaxed:

## Reference
The framework of this project comes from [ahcal-pid](https://github.com/phys-jychen/ahcal-pid). Since their structures are largely different, the definitions of most of the variables have been modified.