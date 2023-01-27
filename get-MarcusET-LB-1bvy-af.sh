echo "
============================================================================================      CYPWare 1.0 script and HTML tool     ==========================================================================================================================
Please cite the use of CYPWare 1.0 and the HTML Tools as follows:
Dixit, V. A.; USN Murty, Bajaj, P.; Blumberger, J.; and Sam P. de Visser. Mechanisms of Electron Transfer Rate Modulations in Cytochrome P450 BM3. J. Phys. Chem. B. 2022 126 (47), 9737-9747. DOI: 10.1021/acs.jpcb.2c03967
https://pubs.acs.org/doi/10.1021/acs.jpcb.2c03967

WELCOME to CYPWare web interface to calculate Marcus ET parameters for Cytochrome P450 BM3 reactions.  
This script requires the following information for successful calculations (also note the dependencies on other software which are required to be in the USER-PATH are given below).
This user-friendly HTML tool (CYPWar-1.0.html) can be used to put all the required information in a text file which CYPWare 1.0 (this script) will read direclty.
If you have these requirements ready, then proceed with the next steps, else perform the docking calculation first and then come back to this CYPWare 1.0 web interface.
1) Name the ligand (without file extension) that was docked into CYP450 active site using Autodock Vina.
2) Name the protein (without file extension) that was used for the Autodock Vina docking simulations.
3) Unix path on your server/workstation where you wish the new files should be created.  You should have write-permissions for that folder. It could be your home directory or any of the sub-directories.
4) Name the folders of the oxidized and reduced states of the Heme center. These directories will be created by the script.
5) Unix path where all the parameter files for the oxidized and reduced states are kept.  These are essential for successfully running the MD simulations.  
	These are kept in a folder called parameters on the GitHub page (in case of query you can post a question to the CYPWare 1.0 development team).
6) Posenumber for the ligand docked pose for which a complex will be created followed by MD simulations and Marcus ET parameter extraction.  
	Currently, CYPWare 1.0 takes Autodock Vina output (pdbqt file) to create the complex, but it can be modified by the user to accept other file formats (which obabel can read).
7) Net molecular charge on the ligand.  This is required to successfully generate ligand parameters with the antechamber program and is required for MD simulations.
The methodology underlying the CYPWare is developed by the PI: Dr. Vaibhav A. Dixit, Asst. Prof., Dept. of Med. Chem., NIPER Guwahati in the Advanced Center of Computer-Aided Drug Design (A-CADD).
The web interface (GUI) is developed in collaboration with C-DAC (Dr. Vinod, Mr. Saurabh, and the team)
CYPWare software, GUI, websites, tools, layouts, and logos are subject to copyright protection and are the exclusive property of the Department of Medicinal Chemistry, NIPER Guwahati. 
The name and logos of the NIPER Guwahati website must not be associated with publicity or business promotion without NIPER G's prior written approval. 
CYPWare 1.0 is available under the creative commons license and is free for academic research groups working in a degree-granting university/institute.  
Any work/report/thesis/research-article/review-article resulting from the use of CYPWare 1.0 should properly cite the software and publication associated with the same.
===========================================================================  Dependencies  ===========================================================================================
CYPWare 1.0 makes use of the following opensource tools, thus these need to be installed first from GitHub, Sourceforge, or as a conda package.
Ensure that dependencies are satisfied before running CYPWare 1.0, else the calculation will not complete as expected.
1) Openbabel 3.1 or higher 
	Openbabel is available as a conda package and can be installed with one of the following commands.
	conda install -c conda-forge openbabel
	conda install -c conda-forge/label/cf202003 openbabel
If you don't have conda, then install it from the main website https://www.anaconda.com/
Instructions for conda installation can be found on its website.
2) AmberTools and Amber18 or higher
	
	AmberTools is a freely available software used for the setup and analysis of MD simulations. It is available from the http://ambermd.org/ website.
	It can also be installed as a conda package with any one of the following command.
	
	conda install -c conda-forge ambertools
	conda install -c conda-forge/label/cf202003 ambertools
	Amber18, 20 or the latest 22 version is a widely used MD engine and includes sander, pmemd, and their serial, parallel (MPI), and GPU (cuda) versions.
	It is available at a reasonable price from Prof. David Case's group at UCSF http://ambermd.org/GetAmber.php#amber
	AMBERHOME directory should be in your path for CYPWare 1.0 to run correctly
3) AmberMdPrep is a wrapper script by Daniel R. Roe used for equilibrating a biomolecular system for MD simulations
	It can be downloaded from https://github.com/drroe/AmberMdPrep
	Untar the AmberMdPrep folder and make sure that the AmberMdPrep.sh script is in your path
	
	AmberMdPrep depends on the latest GitHub version of cpptraj which is available here https://github.com/Amber-MD/cpptraj
	The GitHub version of cpptraj should be installed in a separate directory and sourced before running AmberMdPrep.  CYPWare 1.0 will source cpptraj automatically if the path is set correctly.
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
5) Extracting Marcus ET parameters
	This can be done by calling additional scripts on the Linux terminal.
	bash get-Marcus-ET-parm-LB.sh file (for ligand-bound state) - This script
	bash get-Marcus-ET-parm-LF.sh file (for ligand-free state)
For commercial usage of CYPWare 1.0, please contact the PI at vaibhavadixit@gmail.com or vaibhav@niperguwahati.in
============================================================================================================================================================================================================================================= "

source /home/$USER/.bashrc
source $AMBERHOME/amber.sh

cwd=$(pwd)
export $(xargs <$1)

filecheck () {
	for i in $(ls $dirpath/$1/*fe*.out); do
		if [[ -s $i ]]; 
		then
			echo "output files $i exists thus moving to check other files and extracting data"
		else
			echo "output file $i is missing, check if you ran the MD and SP jos correctly and do the needful"
		fi
	done
}

filecheck $oxddir	
filecheck $reddir

abs() {
  declare -i _value
  _value=$1
  (( _value < 0 )) && _value=$(( _value * -1 ))
  printf "%d\n" $_value
}

> $dirpath/$oxddir/$protein-$ligand-vertEs.txt
> $dirpath/$oxddir/AvgDE_$oxddir.txt
> $dirpath/$oxddir/AvgStdDE_$oxddir.txt
> $dirpath/$oxddir/AvgDEprot_$oxddir.txt
> $dirpath/$oxddir/AvgStdDEprot_$oxddir.txt
> $dirpath/$reddir/$protein-$ligand-vertEs.txt
> $dirpath/$reddir/AvgDE_$reddir.txt
> $dirpath/$reddir/AvgStdDE_$reddir.txt
> $dirpath/$reddir/AvgDEprot_$reddir.txt
> $dirpath/$reddir/AvgStdDEprot_$reddir.txt

calverE () {

        cd $dirpath/$1
	echo "printing sto1D for $1"
        sto1D=$(grep -A1 NSTEP $dirpath/$1/$protein-$ligand-solv-$3.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)
	echo "printing sto2D for $1"
        sto2D=$(grep -A1 NSTEP $dirpath/$1/$protein-$ligand-solv-$4.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)

	echo "printing pto1D for $1"
        pto1D=$(grep -A1 NSTEP $dirpath/$1/$protein-$ligand-prot-$3.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)
	echo "printing pto2D for $1"
        pto2D=$(grep -A1 NSTEP $dirpath/$1/$protein-$ligand-prot-$4.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)

        echo "$sto1D $sto2D $pto1D $pto2D $5" >> $dirpath/$1/$protein-$ligand-vertEs.txt
	echo "printing sto1E for $1"
        sto1E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $2}' | tail -1)
        echo "printing stdev1E  for $1"
        stdev1E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $3}' | tail -1)
        echo "printing sto2E for $1"
        sto2E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $5}' | tail -1)
        echo "printing sto1E for $1"
        stdev2E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $6}' | tail -1)

        echo "printing pto1E for $1"
	pto1E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $8}' | tail -1)
        echo "printing ptdev1E for $1"
	ptdev1E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $9}' | tail -1)
        echo "printing pto2E for $1"
	pto2E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $11}' | tail -1)
        echo "printing ptdev2E for $1"
	ptdev2E=$(cat $dirpath/$1/$protein-$ligand-vertEs.txt | awk '{print $12}' | tail -1)

        echo "printing AvgDE for $1"
	AvgDE=$(bc -l <<< "scale=6; ( ($sto1E - $sto2E) * 0.043 )" )
	echo "AvgDE $AvgDE $5" >> $dirpath/$1/AvgDE_$1.txt
        echo "printing AvgStdDE for $1"
	AvgStdDE=$(bc -l <<< "scale=6; ( ($stdev1E - $stdev2E) * 0.043 )" )
	echo "AvgStdDE $AvgStdDE $5" >> $dirpath/$1/AvgStdDE_$1.txt

        echo "printing AVgDEprot for $1"
	AvgDEprot=$(bc -l <<< "scale=6; ( ($pto1E - $pto2E) * 0.043 )" )
	echo "AvgDEprot $AvgDEprot $5">> $dirpath/$1/AvgDEprot_$1.txt
        echo "printing AvgStdDEprot for $1"
	AvgStdDEprot=$(bc -l <<< "scale=6; ( ($ptdev1E - $ptdev2E) * 0.043 )" )
	echo "AvgStdDEprot $AvgStdDEprot $5">> $dirpath/$1/AvgStdDEprot_$1.txt

}

calverE $oxddir $reddir fe3crdfe3prm fe3crdfe2prm prm
calverE $reddir $oxddir fe2crdfe2prm fe2crdfe3prm prm

calverE $oxddir $reddir fe3crdfe3prm-rp fe3crdfe2prm-rp -rp
calverE $reddir $oxddir fe2crdfe2prm-rp fe2crdfe3prm-rp -rp

calverE $oxddir $reddir fe3crdfe3prm-rp1 fe3crdfe2prm-rp1 -rp1
calverE $reddir $oxddir fe2crdfe2prm-rp1 fe2crdfe3prm-rp1 -rp1

echo "Printing Marcus ET parameters for $protien $ligand complex " > $dirpath/$ligand-MarcusET-parm.txt
echo "MD-traj AvgDEa AvgStdDEa AvgDEb AvgStdDEb lambdaET lambdaStdET DGET StdDGET lambdaETprot lambdaStdETprot cofacdis cofacdisStd" >> $dirpath/$ligand-MarcusET-parm.txt
calMarcusETparm () {
	echo "printing-for-trajectory $1"
	AvgDEred=$(cat $dirpath/$reddir/AvgDE_$reddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEred=$(cat $dirpath/$reddir/AvgStdDE_$reddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEred=$(echo $AvgStdDEred | awk '{print ($1>=0)? $1:0-$1}'  )
	AvgDEoxd=$(cat $dirpath/$oxddir/AvgDE_$oxddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEoxd=$(cat $dirpath/$oxddir/AvgStdDE_$oxddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEoxd=$(echo $AvgStdDEoxd | awk '{print ($1>=0)? $1:0-$1}'  )
	lambdaET=$(bc -l <<< "scale=6; ((- $AvgDEred - $AvgDEoxd)/2)")
	lambdaETStd=$(bc -l <<< "scale=6; ((- $AvgStdDEred - $AvgStdDEoxd))")
	lambdaETStd=$(echo $lambdaETStd | awk '{print ($1>=0)? $1:0-$1}'  )
	DGET=$(bc -l <<< "scale=6; ((- $AvgDEred + $AvgDEoxd)/2)")
	StdDGET=$(bc -l <<< "scale=6; ((- $AvgStdDEred - $AvgStdDEoxd))")
	StdDGET=$(echo $StdDGET | awk '{print ($1>=0)? $1:0-$1}' )

	AvgDEprotred=$(cat $dirpath/$reddir/AvgDEprot_$reddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEprotred=$(cat $dirpath/$reddir/AvgStdDEprot_$reddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEprotred=$(echo $AvgStdDEprotred | awk '{print ($1>=0)? $1:0-$1}'  )
	AvgDEprotoxd=$(cat $dirpath/$oxddir/AvgDEprot_$oxddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEprotoxd=$(cat $dirpath/$oxddir/AvgStdDEprot_$oxddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgStdDEprotoxd=$(echo $AvgStdDEprotoxd | awk '{print ($1>=0)? $1:0-$1}')
	lambdaETprot=$(bc -l <<< "scale=6; ((- $AvgDEprotred - $AvgDEprotoxd)/2)")
	lambdaStdETprot=$(bc -l <<< "scale=6; ((- $AvgStdDEprotred - $AvgStdDEprotoxd))")
	lambdaStdETprot=$(echo $lambdaStdETprot | awk '{print ($1>=0)? $1:0-$1}'  )
	AvgDEred=$(bc -l <<< "scale=6; (- $AvgDEred)")
#	DGETprot=$(bc -l <<< "scale=6; (($AvgDEprotred + $AvgDEprotoxd)/2)"
	filename=$1
	if [[ $filename == prm ]]; 
	then
	        oxdcofacdis=$(tail -n +2 $dirpath/$oxddir/$protein-$ligand-FMN-HEMcontacts.dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        redcofacdis=$(tail -n +2 $dirpath/$reddir/$protein-$ligand-FMN-HEMcontacts.dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        oxdcofacdisstd=$(tail -n +2 $dirpath/$oxddir/$protein-$ligand-FMN-HEMcontacts.dat | awk '{print $4}' | st | tail -1 | awk '{print $6}')
		redcofacdisstd=$(tail -n +2 $dirpath/$reddir/$protein-$ligand-FMN-HEMcontacts.dat | awk '{print $4}' | st | tail -1 | awk '{print $6}')
		avgcofacdis=$(bc -l <<< "scale=6; ( ( $oxdcofacdis + $redcofacdis )/2)" | xargs printf "%6.6f")
		avgcofacdisstd=$(bc -l <<< "scale=6; ( ( $oxdcofacdisstd + $redcofacdisstd ))")
		echo "$1 $AvgDEoxd $AvgStdDEoxd $AvgDEred $AvgStdDEred $lambdaET $lambdaETStd $DGET $StdDGET $lambdaETprot $lambdaStdETprot $avgcofacdis $avgcofacdisstd" >> $dirpath/$ligand-MarcusET-parm.txt
	else
		oxdcofacdis=$(tail -n +2 $dirpath/$oxddir/$protein-$ligand-FMN-HEMcontacts"$1".dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        redcofacdis=$(tail -n +2 $dirpath/$reddir/$protein-$ligand-FMN-HEMcontacts"$1".dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
                oxdcofacdisstd=$(tail -n +2 $dirpath/$oxddir/$protein-$ligand-FMN-HEMcontacts"$1".dat | awk '{print $4}' | st | tail -1 | awk '{print $6}')
                redcofacdisstd=$(tail -n +2 $dirpath/$reddir/$protein-$ligand-FMN-HEMcontacts"$1".dat | awk '{print $4}' | st | tail -1 | awk '{print $6}')
	        avgcofacdis=$(bc -l <<< "scale=6; ( ( $oxdcofacdis + $redcofacdis )/2)" | xargs printf "%6.6f")
                avgcofacdisstd=$(bc -l <<< "scale=6; ( ( $oxdcofacdisstd + $redcofacdisstd ))")
		echo "$1 $AvgDEoxd $AvgStdDEoxd $AvgDEred $AvgStdDEred $lambdaET $lambdaETStd $DGET $StdDGET $lambdaETprot $lambdaStdETprot $avgcofacdis $avgcofacdisstd" >> $dirpath/$ligand-MarcusET-parm.txt
	fi
	
}

calMarcusETparm prm
calMarcusETparm -rp
calMarcusETparm -rp1
more $dirpath/$ligand-MarcusET-parm.txt

