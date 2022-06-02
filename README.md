# CYPWare-1.0
This repository contains a user friendly web-tool (cypwarecode.html) to bed used with CYPWare-1.0 program for calculating CYP450 ET parameters
**-- CYPWare 1.0 script and HTML Tool --**

WELCOME to CYPWare web interface to calculate Marcus ET parameters for Cytochrome P450 BM3 reactions.

This script requires the following information for successful calculations (also note the dependencies on other software which are required to be in the USER-PATH are given below).
This user-friendly HTML tool (cypwarecode.html) can be used to put all the required information in a text file which CYPWare 1.0 (this script) will read direclty.

If you have these requirements ready, then proceed with the next steps, else perform the docking calculation first and then come back to this CYPWare 1.0 web interface.

1) Name the ligand (without file extension) that was docked into CYP450 active site using Autodock Vina.
2) Name the protein (without file extension) that was used for the Autodock Vina docking simulations.
3) Unix path on your server/workstation where you wish the new files should be created. You should have write-permissions for that folder. It could be your home directory or any of the sub-directories.
4) Name the folders of the oxidized and reduced states of the Heme center. These directories will be created by the script.
5) Unix path where all the parameter files for the oxidized and reduced states are kept (the two sub-folders within the BM3-parameters folder). These are essential for successfully running the MD simulations. These are kept in a folder called BM3-parameters on this GitHub page (in case of query you can post a question to the CYPWare 1.0 development team).
6) Posenumber for the ligand docked pose for which a complex will be created followed by MD simulations and Marcus ET parameter extraction. Currently, CYPWare 1.0 takes Autodock Vina output (pdbqt file) to create the complex, but it can be modified by the user to accept other file formats (which obabel can read).
7) Net molecular charge on the ligand. This is required to successfully generate ligand parameters with the antechamber program and is required for MD simulations.

The methodology underlying the CYPWare is developed by the PI: Dr. Vaibhav A. Dixit, Asst. Prof., Dept. of Med. Chem., NIPER Guwahati in the Advanced Center of Computer-Aided Drug Design (A-CADD). The web interface (GUI) is developed in collaboration with C-DAC (Mr. Saurabh G. , Mr. Nikhil R. , Dr. Rajendra J.)
CYPWare software, GUI, websites, tools, layouts, and logos are subject to copyright protection and are the exclusive property of the Department of Medicinal Chemistry, NIPER Guwahati. The name and logos of the NIPER Guwahati website must not be associated with publicity or business promotion without NIPER G's prior written approval.

CYPWare 1.0 is available under the creative commons license and is free for academic research groups working in a degree-granting university/institute.
Any work/report/thesis/research-article/review-article resulting from the use of CYPWare 1.0 should properly cite the software and publication associated with the same.

**-- Dependencies --**

CYPWare 1.0 makes use of the following opensource tools, thus these need to be installed first from GitHub, Sourceforge, or as a conda package.
Ensure that dependencies are satisfied before running CYPWare 1.0, else the calculation will not complete as expected.

1) Openbabel 3.1 or higher

Openbabel is available as a conda package and can be installed with one of the following commands.

conda install -c conda-forge openbabel

conda install -c conda-forge/label/cf202003 openbabel

If you don't have conda, then install it from the main website https://www.anaconda.com/
Instructions for conda installation can be found on its website.

2) AmberTools and Amber18 or higher


AmberTools is a freely available software used for the setup and analysis of MD simulations. It is available from the http://ambermd.org/ website. It can also be installed as a conda package with any one of the following command.
conda install -c conda-forge ambertools

conda install -c conda-forge/label/cf202003 ambertools

Amber18, 20 or the latest 22 version is a widely used MD engine and includes sander, pmemd, and their serial, parallel (MPI), and GPU (cuda) versions.
It is available at a reasonable price from Prof. David Case's group at UCSF http://ambermd.org/GetAmber.php#amber

AMBERHOME directory should be in your path for CYPWare 1.0 to run correctly

3) AmberMdPrep is a wrapper script by Daniel R. Roe used for equilibrating a biomolecular system for MD simulations

It can be downloaded from https://github.com/drroe/AmberMdPrep
Untar the AmberMdPrep folder and make sure that the AmberMdPrep.sh script is in your path

AmberMdPrep depends on the latest GitHub version of cpptraj which is available here https://github.com/Amber-MD/cpptraj
The GitHub version of cpptraj should be installed in a separate directory and sourced before running AmberMdPrep. CYPWare 1.0 will source cpptraj automatically if the path is set correctly.

4) Statistical tool st
A simple statistical tool available on GitHub is used to extract the averages and standard deviations in vertical energy gaps.
This is available at https://github.com/nferraz/st
st can be installed using the following commands

git clone https://github.com/nferraz/st

cd perl5

perl Makefile.PL

make

make test

make install

5) Running CYPWare web tool followed by CYPWare 1.0 program (script)
a) First download the CYPWare package form this repository and extract the contents into a folder of your choice.
b) Then double-click (start) the HTML tool (cypwarecode.html) to open a browser window.
c) Enter ligand, protein file names (without extension) and give directory path and variable values as instructed.
d) Save the file in the location where you want to run the script.  A text file will be saved.
e) Copy the CYPWare 1.0 program (script) into a folder where ligand and protein pdbqt files along with the text file created in step d above.
f) Then run the CYPWare 1.0 program using the following bash commands for liand-bound (LB) or ligand-free (LF) CYP450 BM3 states.

bash CYPWare-1.0-LB-1bvy-af-Marcus-ET-MD.sh file (for ligand-bound state, use file saved in step d)

bash CYPWare-1.0-LF-1bvy-af-Marcus-ET-MD.sh file (for ligand-free state, use file saved in step d)

6) Extracting Marcus ET parameters
This can be done by calling additional scripts on the Linux terminal.

bash get-Marcus-ET-parm-LB.sh file (for ligand-bound state)

bash get-Marcus-ET-parm-LF.sh file (for ligand-free state)

For commercial usage of CYPWare 1.0, please contact the PI at vaibhavadixit@gmail.com or vaibhav@niperguwahati.in

**Funding Sources**
VAD acknowledges the financial support from the National Supercomputing Mission (NSM), Department of Science and Technology (DST), New Delhi 
