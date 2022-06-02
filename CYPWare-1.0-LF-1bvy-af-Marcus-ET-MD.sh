module load cuda/10.1
module load new_amber/20
source /home/$USER/.bashrc
source $AMBERHOME/amber.sh
module load openbabel/2.4.0
cwd=$(pwd)
export $(xargs <$1)
echo "Give a unique name to this run to distinguish it from other Ligand-Free protein MD runs. 
Consider what difference you want to introduce in this run e.g. different salt concentration, type etc. 
Make appropriate changes in the tleap input file below and rerun the script.  
If its not different then seriously consider why you are running it again"

mkdir $dirpath/$reddir
mkdir $dirpath/$oxddir


#echo "Give the name of the ligand without the file extension to be used for complex formation. Be sure that this ligand structure file came from vina docking in the same protein"
#read ligand
echo "Give the name of the protein without the file extension to be used for complex formation. Be sure that this protein structure file was used for vina docking"
echo "Changing directory to Ferrous or the reduced state i.e. $dirpath/$reddir"
cp $protein.pdb $dirpath/$reddir
cp $protein.pdb $dirpath/$oxddir

############################################################
#echo "parameterizing ligand and generating the protein complex"
#babel $ligand.pdbqt -O $ligand.pdb -h -m 
#rm $ligand2.pdb $ligand3.pdb $ligand4.pdb $ligand5.pdb $ligand6.pdb $ligand7.pdb $ligand8.pdb $ligand9.pdb
#grep ATOM $ligand1.pdb > $ligand.pdb
#echo "Give the expeced net-charge on the molecule as an integer e.g. -1, 0 or 1"
#babel test-$ligand.pdb -O test.mol2 --partialcharge gasteiger
#molchrg=$(egrep -A1000 ATOM test.mol2 | grep -v ROOT | grep LIG | sed -e "1d" | awk '{print $9}' | st | awk '{print $4}' | tail -1 | xargs printf "%1.0f")
#read molchrg
#$AMBERHOME/bin/antechamber -i $ligand.pdb -fi pdb -o $ligand.mol2 -fo mol2 -c bcc -s 2 -nc $molchrg
#$AMBERHOME/bin/parmchk2 -i $ligand.mol2 -f mol2 -o $ligand.frcmod
#$AMBERHOME/bin/antechamber -i $ligand.mol2 -fi mol2 -o $ligand.ac -fo ac 
#$AMBERHOME/bin/antechamber -i $ligand.ac -fi ac -o $ligand-test.pdb -fo pdb 
#cp test-$ligand.pdb test.pdb
#sed 's/[[:space:]]*[^[:space:]]*$//' $ligand-test.pdb | sed 's/[[:space:]]*[^[:space:]]*$//'  > $ligand.pdb
#rm test.pdb
#cat $protein.pdb $ligand.pdb > $protein-$ligand.pdb
#pdb4amber -i $protein-$ligand.pdb -o $protein-$ligand-1.pdb
############################################################

echo "Writing tleap input file for the ferric or oxidized state"
echo '
source leaprc.protein.ff14SB
source leaprc.gaff2
addAtomTypes {' > $dirpath/$oxddir/$protein-tleap.in
echo '        { "M1"  "Fe" "sp3" }' >> $dirpath/$oxddir/$protein-tleap.in
echo '        { "Y1"  "S" "sp3" }' >> $dirpath/$oxddir/$protein-tleap.in
echo '        { "Y2"  "N" "sp3" }' >> /$dirpath/$oxddir/$protein-tleap.in
echo '        { "Y3"  "N" "sp3" }' >> $dirpath/$oxddir/$protein-tleap.in
echo '        { "Y4"  "N" "sp3" }' >> $dirpath/$oxddir/$protein-tleap.in
echo '        { "Y5"  "N" "sp3" }' >> $dirpath/$oxddir/$protein-tleap.in
echo '}' >> $dirpath/$oxddir/$protein-tleap.in
echo "CM1 = loadmol2 $oxdparm/CM1.mol2
HM1 = loadmol2 $oxdparm/HM1.mol2
FE1 = loadmol2 $oxdparm/FE1.mol2
FMN = loadmol2 $oxdparm/FMN-sq-H.mol2
#LIG = loadmol2 $dirpath/$oxddir/$ligand.mol2
loadamberparams $oxdparm/HEM.frcmod
#loadamberparams $dirpath/$oxddir/$ligand.frcmod
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams $oxdparm/4ZF6_mcpbpy.frcmod
loadamberparams $oxdparm/FMN-sq-H.frcmod
source leaprc.water.tip3p
mol = loadpdb $dirpath/$oxddir/$protein.pdb
bond mol.381.SG mol.613.FE
bond mol.612.NA mol.613.FE
bond mol.612.NB mol.613.FE
bond mol.612.NC mol.613.FE
bond mol.612.ND mol.613.FE
bond mol.380.C mol.381.N
bond mol.381.C mol.382.N
savepdb mol $dirpath/$oxddir/$protein-dry.pdb
saveamberparm mol $dirpath/$oxddir/$protein-dry.prmtop $dirpath/$oxddir/$protein-dry.inpcrd
charge mol
solvatebox mol TIP3PBOX 10.0
addions mol K+ 145 
addions mol Cl- 116
charge mol
savepdb mol $dirpath/$oxddir/$protein-solv.pdb
saveamberparm mol $dirpath/$oxddir/$protein-solv.prmtop $dirpath/$oxddir/$protein-solv.inpcrd
quit" >> $dirpath/$oxddir/$protein-tleap.in
echo "tleap input file has been created for the ferric or oxidized state. Please check if the system charged is 0, if not make appropriate changes in the type and number of ions to be added in the addions file and rerun this main script"

echo "Writing tleap input file for the ferrous or reduced state"
echo '
source leaprc.protein.ff14SB
source leaprc.gaff2
addAtomTypes {' > $dirpath/$reddir/$protein-tleap.in
echo '        { "M1"  "Fe" "sp3" }' >> $dirpath/$reddir/$protein-tleap.in
echo '        { "Y1"  "S" "sp3" }' >> $dirpath/$reddir/$protein-tleap.in
echo '        { "Y2"  "N" "sp3" }' >> /$dirpath/$reddir/$protein-tleap.in
echo '        { "Y3"  "N" "sp3" }' >> $dirpath/$reddir/$protein-tleap.in
echo '        { "Y4"  "N" "sp3" }' >> $dirpath/$reddir/$protein-tleap.in
echo '        { "Y5"  "N" "sp3" }' >> $dirpath/$reddir/$protein-tleap.in
echo '}' >> $dirpath/$reddir/$protein-tleap.in
echo "CM1 = loadmol2 $redparm/CM1.mol2
HM1 = loadmol2 $redparm/HM1.mol2
FE1 = loadmol2 $redparm/FE1.mol2
FMN = loadmol2 $redparm/FMN-ox-H.mol2
#LIG = loadmol2 $dirpath/$reddir/$ligand.mol2
loadamberparams $redparm/HEM.frcmod
#loadamberparams $dirpath/$reddir/$ligand.frcmod
loadamberparams frcmod.ions1lm_126_tip3p
loadamberparams $redparm/4ZF6_mcpbpy.frcmod
loadamberparams $redparm/FMN-ox-H.frcmod
source leaprc.water.tip3p
mol = loadpdb $dirpath/$reddir/$protein.pdb
bond mol.381.SG mol.613.FE
bond mol.612.NA mol.613.FE
bond mol.612.NB mol.613.FE
bond mol.612.NC mol.613.FE
bond mol.612.ND mol.613.FE
bond mol.380.C mol.381.N
bond mol.381.C mol.382.N
savepdb mol $dirpath/$reddir/$protein-dry.pdb
saveamberparm mol $dirpath/$reddir/$protein-dry.prmtop $dirpath/$reddir/$protein-dry.inpcrd
charge mol
solvatebox mol TIP3PBOX 10.0
addions mol K+ 145
addions mol Cl- 116
charge mol
savepdb mol $dirpath/$reddir/$protein-solv.pdb
saveamberparm mol $dirpath/$reddir/$protein-solv.prmtop $dirpath/$reddir/$protein-solv.inpcrd
quit" >> $dirpath/$reddir/$protein-tleap.in
echo "tleap input file has been created for the ferric or oxidized state. Please check if the system charged is 0 if not make appropriate changes in the type and number of ions to be added in the addions file and rerun this main script"

echo "Writing cpptraj input files"

writemdfiles () {
	echo "
parm $protein-solv.prmtop
trajin $protein-solv-prod.nc 
autoimage
center
strip :WAT,Cl-,Na+,K+
trajout $protein-prot.nc 
trajout $protein-prot.rst7 onlyframes 1 
" > $dirpath/$1/protein-reorgE-cpptraj.in 

	echo "
parm $protein-solv.prmtop
trajin $protein-solv-prod-rp.nc 
autoimage
center
strip :WAT,Cl-,Na+,K+
trajout $protein-prot-rp.nc 
trajout $protein-prot-rp.rst7 onlyframes 1
" > $dirpath/$1/protein-reorgE-cpptraj1.in

	echo "
parm $protein-solv.prmtop
trajin $protein-solv-prod-rp1.nc 
autoimage
center
strip :WAT,Cl-,Na+,K+
trajout $protein-prot-rp1.nc 
trajout $protein-prot-rp1.rst7 onlyframes 1
" > $dirpath/$1/protein-reorgE-cpptraj2.in

	echo "Writing sp-MD input file"
	echo "single point energy calculations of ferric state coordinates using ferrous parameters
&cntrl
  ntx=1, ntpr=1, ntwx=1,
  imin=5, maxcyc=1, !Single-point energy calculation on each frame
  ntb=1, ntp=0, ntc=2, !periodic
  cut=9, !Calculate all solute-solute interactions
/
" > $dirpath/$1/MD-sp-periodic.mdin

	echo "Writing MD input file for a 20 ns MD production run"
	echo "Production Stage in Explicit Solvent 20 ns
&cntrl
  ntt=3,           ! Temperature scaling (=3, Langevin dynamics)
  gamma_ln=2.0,    ! Collision frequency of the Langevin dynamics in ps-1
  ntc=2,           ! SHAKE constraints (=2, hydrogen bond lengths constrained)
  ntf=2,           ! Force evaluation (=2, hydrogen bond interactions omitted)
  ntb=1,           ! Boundaries (=1, constant volume)
  cut=10.0,        ! Cutoff
  dt=0.002,        ! The time step in picoseconds
  nstlim=10000000, ! Number of MD steps to be performed
  ig=-1,           ! Random seed (=-1, get a number from current date and time)
  ntwr=5000,      ! Restart file written every ntwr steps
  ntwx=5000,      ! Trajectory file written every ntwx steps
  ntpr=5000,      ! The mdout and mdinfo files written every ntpr steps
  ioutfm=1,        ! Trajectory file format (=1, Binary NetCDF)
  iwrap=1,         ! Translate water molecules into the original simulation box
  igb=0,           ! GB model (=0, explicit solvent)
  irest=1,         ! Flag to restart the simulation
  ntx=5,           ! Initial condition (=5, coord. and veloc. read from the inpcrd file)
/
" > $dirpath/$1/prod20ns.in

#Ferrousdir=$(echo Ferrous-$protein-$ligand | sed 's/_out//g')

	echo "Writing sp-MD input files for the ferric/ferrous or oxidized/reduced state to be submitted by task spoolar ts"
	echo "
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-solv-$3.out -c $protein-solv-prod.rst7 -ref $protein-solv-prod.rst7 -r restrt -y $protein-solv-prod.nc -p $protein-solv.prmtop 
wait 
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-solv-$4.out -c $protein-solv-prod.rst7 -ref $protein-solv-prod.rst7 -r restrt -y $protein-solv-prod.nc -p ../$2/$protein-solv.prmtop " > $dirpath/$1/sp-MD.sh

	echo "
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-solv-$3-rp.out -c $protein-solv-prod-rp.rst7 -ref $protein-solv-prod-rp.rst7 -r restrt -y $protein-solv-prod-rp.nc -p $protein-solv.prmtop 
wait 
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-solv-$4-rp.out -c $protein-solv-prod-rp.rst7 -ref $protein-solv-prod-rp.rst7 -r restrt -y $protein-solv-prod-rp.nc -p ../$2/$protein-solv.prmtop " > $dirpath/$1/sp-MD1.sh

	echo "
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-solv-$3-rp1.out -c $protein-solv-prod-rp1.rst7 -ref $protein-solv-prod-rp1.rst7 -r restrt -y $protein-solv-prod-rp1.nc -p $protein-solv.prmtop 
wait 
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-solv-$4-rp1.out -c $protein-solv-prod-rp1.rst7 -ref $protein-solv-prod-rp1.rst7 -r restrt -y $protein-solv-prod-rp1.nc -p ../$2/$protein-solv.prmtop " > $dirpath/$1/sp-MD2.sh

	echo "
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-prot-$3.out -c $protein-prot.rst7 -ref $protein-prot.rst7 -r restrt -y $protein-prot.nc -p $protein-dry.prmtop 
wait 
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-prot-$4.out -c $protein-prot.rst7 -ref $protein-prot.rst7 -r restrt -y $protein-prot.nc -p ../$2/$protein-dry.prmtop " > $dirpath/$1/sp-MD-prot.sh

	echo "
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-prot-$3-rp.out -c $protein-prot-rp.rst7 -ref $protein-prot-rp.rst7 -r restrt -y $protein-prot-rp.nc -p $protein-dry.prmtop 
wait 
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-prot-$4-rp.out -c $protein-prot-rp.rst7 -ref $protein-prot-rp.rst7 -r restrt -y $protein-prot-rp.nc -p ../$2/$protein-dry.prmtop " > $dirpath/$1/sp-MD-prot1.sh

	echo "
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-prot-$3-rp1.out -c $protein-prot-rp1.rst7 -ref $protein-prot-rp1.rst7 -r restrt -y $protein-prot-rp1.nc -p $protein-dry.prmtop 
wait 
mpirun -np 12 sander.MPI -O -i MD-sp-periodic.mdin -o $protein-prot-$4-rp1.out -c $protein-prot-rp1.rst7 -ref $protein-prot-rp1.rst7 -r restrt -y $protein-prot-rp1.nc -p ../$2/$protein-dry.prmtop " > $dirpath/$1/sp-MD-prot2.sh

	echo "
parm $protein-solv.prmtop 
trajin $protein-solv-prod.nc 
autoimage
center
reference final.1.ncrst [GB1]

rms ToGB1 ref [GB1] :1-611@CA,N,C,O :1-611@CA,N,C,O out rmsd-$protein-solv-prod.dat" > $dirpath/$1/rmsd-cpptraj.in

	echo "
parm $protein-solv.prmtop 
trajin $protein-solv-prod-rp.nc 
autoimage
center
reference final.1.ncrst [GB1]

rms ToGB1 ref [GB1] :1-611@CA,N,C,O :1-611@CA,N,C,O out rmsd-$protein-solv-prod-rp.dat" > $dirpath/$1/rmsd-cpptraj1.in

	echo "
parm $protein-solv.prmtop 
trajin $protein-solv-prod-rp1.nc 
autoimage
center
reference final.1.ncrst [GB1]

rms ToGB1 ref [GB1] :1-611@CA,N,C,O :1-611@CA,N,C,O out rmsd-$protein-solv-prod-rp1.dat" > $dirpath/$1/rmsd-cpptraj2.in

	echo "
parm $protein-solv.prmtop
trajin $protein-solv-prod.nc
autoimage 
center
reference final.1.ncrst

nativecontacts name $protein :1-439 :459-611 byresidue out $protein-contacts.dat mindist maxdist distance 7.0 first map mapout $protein-resmap.gnu contactpdb $protein-contactmap.pdb series seriesout $protein-nativecontacts.dat writecontacts $protein-writecontacts.dat
nativecontacts name $protein-HEM :612&!@H= :459-611&!@H= byresidue out $protein-HEMcontacts.dat mindist maxdist distance 3.0 first map mapout $protein-HEMresmap.gnu contactpdb $protein-HEMcontactmap.pdb series seriesout $protein-HEMnativecontacts.dat writecontacts $protein-HEMwritecontacts.dat
nativecontacts name $protein-FMN :614&!@H= :1-439&!@H= byresidue out $protein-FMNcontacts.dat mindist maxdist distance 3.0 first map mapout $protein-FMNresmap.gnu contactpdb $protein-FMNcontactmap.pdb series seriesout $protein-FMNnativecontacts.dat writecontacts $protein-FMNwritecontacts.dat
nativecontacts name $protein-FMN-HEM :612&!@H= :614&!@H= byresidue out $protein-FMN-HEMcontacts.dat mindist maxdist distance 3.0 first map mapout $protein-FMN-HEMresmap.gnu contactpdb $protein-FMN-HEMcontactmap.pdb series seriesout $protein-FMN-HEMnativecontacts.dat writecontacts $protein-FMN-HEMwritecontacts.dat " > $dirpath/$1/contacts.in

	echo "
parm $protein-solv.prmtop
trajin $protein-solv-prod-rp.nc
autoimage 
center
reference final.1.ncrst

nativecontacts name $protein :1-439 :459-611 byresidue out $protein-contacts-rp.dat mindist maxdist distance 7.0 first map mapout $protein-resmap-rp.gnu contactpdb $protein-contactmap-rp.pdb series seriesout $protein-nativecontacts-rp.dat writecontacts $protein-writecontacts-rp.dat
nativecontacts name $protein-HEM :612&!@H= :459-611&!@H= byresidue out $protein-HEMcontacts-rp.dat mindist maxdist distance 3.0 first map mapout $protein-HEMresmap-rp.gnu contactpdb $protein-HEMcontactmap-rp.pdb series seriesout $protein-HEMnativecontacts-rp.dat writecontacts $protein-HEMwritecontacts-rp.dat
nativecontacts name $protein-FMN :614&!@H= :1-439&!@H= byresidue out $protein-FMNcontacts-rp.dat mindist maxdist distance 3.0 first map mapout $protein-FMNresmap-rp.gnu contactpdb $protein-FMNcontactmap-rp.pdb series seriesout $protein-FMNnativecontacts-rp.dat writecontacts $protein-FMNwritecontacts-rp.dat
nativecontacts name $protein-FMN-HEM :612&!@H= :614&!@H= byresidue out $protein-FMN-HEMcontacts-rp.dat mindist maxdist distance 3.0 first map mapout $protein-FMN-HEMresmap-rp.gnu contactpdb $protein-FMN-HEMcontactmap-rp.pdb series seriesout $protein-FMN-HEMnativecontacts-rp.dat writecontacts $protein-FMN-HEMwritecontacts-rp.dat " > $dirpath/$1/contacts-rp.in

	echo "
parm $protein-solv.prmtop
trajin $protein-solv-prod-rp1.nc
autoimage 
center
reference final.1.ncrst

nativecontacts name $protein :1-439 :459-611 byresidue out $protein-contacts-rp1.dat mindist maxdist distance 7.0 first map mapout $protein-resmap-rp1.gnu contactpdb $protein-contactmap-rp1.pdb series seriesout $protein-nativecontacts-rp1.dat writecontacts $protein-writecontacts-rp1.dat
nativecontacts name $protein-HEM :612&!@H= :459-611&!@H= byresidue out $protein-HEMcontacts-rp1.dat mindist maxdist distance 3.0 first map mapout $protein-HEMresmap-rp1.gnu contactpdb $protein-HEMcontactmap-rp1.pdb series seriesout $protein-HEMnativecontacts-rp1.dat writecontacts $protein-HEMwritecontacts-rp1.dat
nativecontacts name $protein-FMN :614&!@H= :1-439&!@H= byresidue out $protein-FMNcontacts-rp1.dat mindist maxdist distance 3.0 first map mapout $protein-FMNresmap-rp1.gnu contactpdb $protein-FMNcontactmap-rp1.pdb series seriesout $protein-FMNnativecontacts-rp1.dat writecontacts $protein-FMNwritecontacts-rp1.dat
nativecontacts name $protein-FMN-HEM :612&!@H= :614&!@H= byresidue out $protein-FMN-HEMcontacts-rp1.dat mindist maxdist distance 3.0 first map mapout $protein-FMN-HEMresmap-rp1.gnu contactpdb $protein-FMN-HEMcontactmap-rp1.pdb series seriesout $protein-FMN-HEMnativecontacts-rp1.dat writecontacts $protein-FMN-HEMwritecontacts-rp1.dat " > $dirpath/$1/contacts-rp1.in

}

writemdfiles $oxddir $reddir fe3crdfe3prm fe3crdfe2prm
writemdfiles $reddir $oxddir fe2crdfe2prm fe2crdfe3prm


#echo "give leap input and output file names one after the other and press enter before entering the next file name"
leapinput=$protein-tleap.in
leapoutput=$protein-tleap.out
#echo "give base file names to be used for MD prmtop, inpcrd, out, rst7, mdinfo and trajectory (nc) files"
basefile=$protein-solv
echo $leapinput $leapoutput $basefile
echo "creating MD prmtop inpcrd files for the Ferric or the oxidized state in $dirpath/$oxddir "

echo "#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --time=95:50:20
#SBATCH --job-name=$protein-$ligand-solv
#SBATCH --error=job.%J.err_node_40
#SBATCH --output=job.%J.out_node_40
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
export I_MPI_FABRICS=shm:dapl
hostname

" > gpu-job-head.txt


echo "creating cpu-job head " | tee -a $protein-$ligand-MD-logfile.txt
echo "#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH --time=95:50:20
#SBATCH --job-name=$protein-$ligand-solv
#SBATCH --error=job.%J.err_node_40
#SBATCH --output=job.%J.out_node_40
#SBATCH --partition=cpu
export I_MPI_FABRICS=shm:dapl
hostname

" > cpu-job-head.txt
echo "creating tleap input and sbatch file for $oxddir state" | tee -a $protein-$ligand-MD-logfile.txt
cat cpu-job-head.txt > $dirpath/$oxddir/$protein-$ligand-tleap.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/$protein-$ligand-tleap.batch
echo "/home/apps/amber/tmp/install_gnu/bin/tleap -s -f $dirpath/$oxddir/$leapinput > $dirpath/$oxddir/$leapoutput " >> $dirpath/$oxddir/$protein-$ligand-tleap.batch


oxdtleapjobid=$(sbatch $dirpath/$oxddir/$protein-$ligand-tleap.batch | awk '{print $4}')
echo "sleeping for 2m"
sleep 2m
oxdtleapjobstatus=$(sacct | grep $oxdtleapjobid | awk '{print $6}' | head -1)

echo "creating tleap input and sbatch file for $reddir state" | tee -a $protein-$ligand-MD-logfile.txt
cat cpu-job-head.txt > $dirpath/$reddir/$protein-$ligand-tleap.batch
echo "creating MD prmtop inpcrd files for the Ferrous or the reduced state in $dirpath/$reddir "
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$reddir/$protein-$ligand-tleap.batch
echo "/home/apps/amber/tmp/install_gnu/bin/tleap -s -f $dirpath/$reddir/$leapinput > $dirpath/$reddir/$leapoutput " >> $dirpath/$reddir/$protein-$ligand-tleap.batch

redtleapjobid=$(sbatch $dirpath/$reddir/$protein-$ligand-tleap.batch | awk '{print $4}')
echo "sleeping for 2m"
sleep 2m
redtleapjobstatus=$(sacct | grep $redtleapjobid | awk '{print $6}' | head -1)




echo "creating input and submitting equil jobs" | tee -a $dirpath/$protein-MD-logfile.txt

while [ $redtleapjobstatus == RUNNING ];
do
	sleep 1m
done
echo "finished equiliration of $oxddir structure" | tee -a $dirpath/$protein-MD-logfile.txt

while [ $oxdtleapjobstatus == RUNNING ];
do
	sleep 1m
done
echo "finished equiliration of $reddir structure" | tee -a $dirpath/$protein-MD-logfile.txt
#ts -S 1 
chmod a+x $dirpath/$oxddir/sp-MD*.sh
chmod a+x $dirpath/$reddir/sp-MD*.sh

source /home/vaibhavdixit.niperg/cpptraj-master/cpptraj.sh
echo "changing directory to $oxddir" | tee -a $dirpath/$protein-MD-logfile.txt
cd $dirpath/$oxddir
cp ../gpu-job-head.txt .
echo "creating MDequil sbatch files and submitting the job $oxddir state" | tee -a $dirpath/$protein-MD-logfile.txt
cat gpu-job-head.txt > $dirpath/$oxddir/$protein-MDequil.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/$protein-MDequil.batch
echo "source /home/vaibhavdixit.niperg/cpptraj-master/cpptraj.sh" >> $dirpath/$oxddir/$protein-MDequil.batch
echo "/home/vaibhavdixit.niperg/AmberMdPrep-master/AmberMdPrep.sh -O -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/$basefile.inpcrd --temp 300 --ares HD1 --ares HM1 --ares FE1 --ares CM1 --ares FMN " >> $dirpath/$oxddir/$protein-MDequil.batch
oxdequiljobid=$(sbatch $dirpath/$oxddir/$protein-MDequil.batch | awk '{print $4}')
oxdequiljobstatus=$(sacct | grep $oxdequiljobid | awk '{print $6}' | head -1)

source /home/vaibhavdixit.niperg/cpptraj-master/cpptraj.sh
cd $dirpath/$reddir
cp ../gpu-job-head.txt .
echo "creating MDequil sbatch files and submitting the job $reddir state" | tee -a $dirpath/$protein-MD-logfile.txt
cat gpu-job-head.txt > $dirpath/$reddir/$protein-MDequil.batch
echo "source /home/vaibhavdixit.niperg/cpptraj-master/cpptraj.sh" >> $dirpath/$reddir/$protein-MDequil.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$reddir/$protein-MDequil.batch
echo "/home/vaibhavdixit.niperg/AmberMdPrep-master/AmberMdPrep.sh -O -p $dirpath/$reddir/$basefile.prmtop -c $dirpath/$reddir/$basefile.inpcrd --temp 300 --ares HD1 --ares HM1 --ares FE1 --ares CM1 --ares FMN " >> $dirpath/$reddir/$protein-MDequil.batch
redequiljobid=$(sbatch $dirpath/$reddir/$protein-MDequil.batch | awk '{print $4}')
redequiljobstatus=$(sacct | grep $redequiljobid | awk '{print $6}' | head -1)

echo "sleeping for the equiljobs to finish"

while [[ $(sacct | grep $oxdequiljobid | awk '{print $6}' | head -1 ) == RUNNING ]];
do
	sleep 1m
done
echo "finished equiliration of $oxddir structure" | tee -a $dirpath/$protein-MD-logfile.txt

while [[ $(sacct | grep $redequiljobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 1m
done


#if [[ $oxdMDstatus == "Equilibration success" ]]
#	then
echo "changing directory to $oxddir and submitting 20ns MD jobs in triplicate (prod, prod-rp and prod-rp1) after equilibriation is success" | tee -a $dirpath/$protein-MD-logfile.txt
cd $dirpath/$oxddir
cat gpu-job-head.txt > $dirpath/$oxddir/oxdMD-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/oxdMD-rp-rp1.batch
echo "pmemd.cuda -O -i $dirpath/$oxddir/prod20ns.in -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/final.1.ncrst -ref $dirpath/$oxddir/final.1.ncrst -r $dirpath/$oxddir/$basefile-prod.rst7 -inf $dirpath/$oxddir/$basefile-prod.mdinfo -o $dirpath/$oxddir/$basefile-prod.mdout -x $dirpath/$oxddir/$basefile-prod.nc 
pmemd.cuda -O -i $dirpath/$oxddir/prod20ns.in -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/final.1.ncrst -ref $dirpath/$oxddir/final.1.ncrst -r $dirpath/$oxddir/$basefile-prod-rp.rst7 -inf $dirpath/$oxddir/$basefile-prod-rp.mdinfo -o $dirpath/$oxddir/$basefile-prod-rp.mdout -x $dirpath/$oxddir/$basefile-prod-rp.nc 
pmemd.cuda -O -i $dirpath/$oxddir/prod20ns.in -p $dirpath/$oxddir/$basefile.prmtop -c $dirpath/$oxddir/final.1.ncrst -ref $dirpath/$oxddir/final.1.ncrst -r $dirpath/$oxddir/$basefile-prod-rp1.rst7 -inf $dirpath/$oxddir/$basefile-prod-rp1.mdinfo -o $dirpath/$oxddir/$basefile-prod-rp1.mdout -x $dirpath/$oxddir/$basefile-prod-rp1.nc 
cpptraj.cuda -i $dirpath/$oxddir/contacts.in
cpptraj.cuda -i $dirpath/$oxddir/contacts-rp.in
cpptraj.cuda -i $dirpath/$oxddir/contacts-rp1.in 

" >> $dirpath/$oxddir/oxdMD-rp-rp1.batch 
oxdMDjobid=$(sbatch --dependency=afterok:$oxdequiljobid $dirpath/$oxddir/oxdMD-rp-rp1.batch | awk '{print $4}')
oxdMDjobstatus=$(sacct | grep $oxdMDjobid | awk '{print $6}'  | head -1)


echo "changing directory to $reddir and submitting 20ns MD jobs in triplicate (prod, prod-rp and prod-rp1) after equilibriation is success" | tee -a $dirpath/$protein-MD-logfile.txt
cd $dirpath/$reddir
cat gpu-job-head.txt > $dirpath/$reddir/redMD-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$reddir/redMD-rp-rp1.batch
echo "pmemd.cuda -O -i $dirpath/$reddir/prod20ns.in -p $dirpath/$reddir/$basefile.prmtop -c $dirpath/$reddir/final.1.ncrst -ref $dirpath/$reddir/final.1.ncrst -r $dirpath/$reddir/$basefile-prod.rst7 -inf $dirpath/$reddir/$basefile-prod.mdinfo -o $dirpath/$reddir/$basefile-prod.mdout -x $dirpath/$reddir/$basefile-prod.nc 
pmemd.cuda -O -i $dirpath/$reddir/prod20ns.in -p $dirpath/$reddir/$basefile.prmtop -c $dirpath/$reddir/final.1.ncrst -ref $dirpath/$reddir/final.1.ncrst -r $dirpath/$reddir/$basefile-prod-rp.rst7 -inf $dirpath/$reddir/$basefile-prod-rp.mdinfo -o $dirpath/$reddir/$basefile-prod-rp.mdout -x $dirpath/$reddir/$basefile-prod-rp.nc 
pmemd.cuda -O -i $dirpath/$reddir/prod20ns.in -p $dirpath/$reddir/$basefile.prmtop -c $dirpath/$reddir/final.1.ncrst -ref $dirpath/$reddir/final.1.ncrst -r $dirpath/$reddir/$basefile-prod-rp1.rst7 -inf $dirpath/$reddir/$basefile-prod-rp1.mdinfo -o $dirpath/$reddir/$basefile-prod-rp1.mdout -x $dirpath/$reddir/$basefile-prod-rp1.nc 
cpptraj.cuda -i $dirpath/$reddir/contacts.in
cpptraj.cuda -i $dirpath/$reddir/contacts-rp.in
cpptraj.cuda -i $dirpath/$reddir/contacts-rp1.in
" >> $dirpath/$reddir/redMD-rp-rp1.batch 
redMDjobid=$(sbatch --dependency=afterok:$redequiljobid $dirpath/$reddir/redMD-rp-rp1.batch | awk '{print $4}')
redMDjobstatus=$(sacct | grep $redMDjobid | awk '{print $6}'  | head -1)


while [[ $(sacct | grep $oxdMDjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 5m
done

while [[ $(sacct | grep $redMDjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 5m
done

echo "changing directory to $oxddir and submitting cpptraj, SP MD jobs for MD trajectories (prod, prod-rp and prod-rp1) after production runs is success" | tee -a $dirpath/$protein-MD-logfile.txt
cd $dirpath/$oxddir
cat ../cpu-job-head.txt > $dirpath/$oxddir/cpptraj-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/cpptraj-rp-rp1.batch
echo "cpptraj -i $dirpath/$oxddir/protein-reorgE-cpptraj.in
cpptraj -i $dirpath/$oxddir/protein-reorgE-cpptraj1.in
cpptraj -i $dirpath/$oxddir/protein-reorgE-cpptraj2.in" >> $dirpath/$oxddir/cpptraj-rp-rp1.batch
oxdcpptrajjobid=$(sbatch --dependency=afterok:$oxdMDjobid $dirpath/$oxddir/cpptraj-rp-rp1.batch | awk '{print $4}')


cat ../cpu-job-head.txt > $dirpath/$oxddir/sp-MD-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/sp-MD-rp-rp1.batch
echo "$dirpath/$oxddir/./sp-MD.sh 
$dirpath/$oxddir/./sp-MD1.sh 
$dirpath/$oxddir/./sp-MD2.sh " >> $dirpath/$oxddir/sp-MD-rp-rp1.batch
oxdspMDjobid=$(sbatch --dependency=afterok:$oxdcpptrajjobid $dirpath/$oxddir/sp-MD-rp-rp1.batch | awk '{print $4}')


cat ../cpu-job-head.txt > $dirpath/$oxddir/sp-MD-prot-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$oxddir/sp-MD-prot-rp-rp1.batch
echo "$dirpath/$oxddir/./sp-MD-prot.sh 
$dirpath/$oxddir/./sp-MD-prot1.sh 
$dirpath/$oxddir/./sp-MD-prot2.sh " >> $dirpath/$oxddir/sp-MD-prot-rp-rp1.batch
oxdspMDprotjobid=$(sbatch --dependency=afterok:$oxdcpptrajjobid $dirpath/$oxddir/sp-MD-prot-rp-rp1.batch | awk '{print $4}')



echo "changing directory to $reddir and submitting cpptraj, SP MD jobs for MD trajectories (prod, prod-rp and prod-rp1) after production runs is success" | tee -a $dirpath/$protein-MD-logfile.txt
cd $dirpath/$reddir
cat ../cpu-job-head.txt > $dirpath/$reddir/cpptraj-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$reddir/cpptraj-rp-rp1.batch
echo "cpptraj -i $dirpath/$reddir/protein-reorgE-cpptraj.in
cpptraj -i $dirpath/$reddir/protein-reorgE-cpptraj1.in
cpptraj -i $dirpath/$reddir/protein-reorgE-cpptraj2.in" >> $dirpath/$reddir/cpptraj-rp-rp1.batch
redcpptrajjobid=$(sbatch --dependency=afterok:$redMDjobid $dirpath/$reddir/cpptraj-rp-rp1.batch | awk '{print $4}')


cat ../cpu-job-head.txt > $dirpath/$reddir/sp-MD-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$reddir/sp-MD-rp-rp1.batch
echo "$dirpath/$reddir/./sp-MD.sh 
$dirpath/$reddir/./sp-MD1.sh 
$dirpath/$reddir/./sp-MD2.sh " >> $dirpath/$reddir/sp-MD-rp-rp1.batch
redspMDjobid=$(sbatch --dependency=afterok:$redcpptrajjobid $dirpath/$reddir/sp-MD-rp-rp1.batch | awk '{print $4}')


cat ../cpu-job-head.txt > $dirpath/$reddir/sp-MD-prot-rp-rp1.batch
echo "module load new_amber/20
source $AMBERHOME/amber.sh" >> $dirpath/$reddir/sp-MD-prot-rp-rp1.batch
echo "$dirpath/$reddir/./sp-MD-prot.sh 
$dirpath/$reddir/./sp-MD-prot1.sh 
$dirpath/$reddir/./sp-MD-prot2.sh " >> $dirpath/$reddir/sp-MD-prot-rp-rp1.batch
redspMDprotjobid=$(sbatch --dependency=afterok:$redcpptrajjobid $dirpath/$reddir/sp-MD-prot-rp-rp1.batch | awk '{print $4}')


oxdspMDjobstatus=$(sacct | grep $oxdspMDjobid | awk '{print $6}'  | head -1 )

redspMDjobstatus=$(sacct | grep $redspMDjobid | awk '{print $6}'  | head -1 )

oxdspMDprotjobstatus=$(sacct | grep $oxdspMDprotjobid | awk '{print $6}'  | head -1 )

redspMDprotjobstatus=$(sacct | grep $redspMDprotjobid | awk '{print $6}'  | head -1 )



while [[ $(sacct | grep $oxdspMDjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 5m
done

while [[ $(sacct | grep $redspMDjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 5m
done

while [[ $(sacct | grep $oxdspMDprotjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 2m
done

while [[ $(sacct | grep $redspMDprotjobid | awk '{print $6}'  | head -1 ) == RUNNING ]];
do
	sleep 2m
done



> $dirpath/$oxddir/$protein-$uniquename-vertEs.txt
> $dirpath/$oxddir/AvgDE_$oxddir.txt
> $dirpath/$oxddir/AvgStdDE_$oxddir.txt
> $dirpath/$oxddir/AvgDEprot_$oxddir.txt
> $dirpath/$oxddir/AvgStdDEprot_$oxddir.txt
> $dirpath/$reddir/$protein-vertEs.txt
> $dirpath/$reddir/AvgDE_$reddir.txt
> $dirpath/$reddir/AvgStdDE_$reddir.txt
> $dirpath/$reddir/AvgDEprot_$reddir.txt
> $dirpath/$reddir/AvgStdDEprot_$reddir.txt

calverE () {

        cd $dirpath/$1
        sto1D=$(grep -A1 NSTEP $dirpath/$1/$protein-solv-$3.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)
        sto2D=$(grep -A1 NSTEP $dirpath/$1/$protein-solv-$4.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)

        pto1D=$(grep -A1 NSTEP $dirpath/$1/$protein-prot-$3.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)
        pto2D=$(grep -A1 NSTEP $dirpath/$1/$protein-prot-$4.out | grep -v NSTEP | sed '/--/d' | awk '{print $2}' | st |  awk '{print $1, $5, $6}' | tail -n+2 | head -1)

        echo "$sto1D $sto2D $pto1D $pto2D $5" >> $dirpath/$1/$protein-vertEs.txt
        sto1E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $2}' | tail -1)
        stdev1E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $3}' | tail -1)
        sto2E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $5}' | tail -1)
        stdev2E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $6}' | tail -1)

	pto1E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $8}' | tail -1)
	ptdev1E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $9}' | tail -1)
	pto2E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $11}' | tail -1)
	ptdev2E=$(cat $dirpath/$1/$protein-vertEs.txt | awk '{print $12}' | tail -1)

	AvgDE=$(bc -l <<< "scale=6; ( ($sto1E - $sto2E) * 0.043 )" | xargs printf "%6.6f")
	echo "AvgDE $AvgDE $5" >> $dirpath/$1/AvgDE_$1.txt
	AvgStdDE=$(bc -l <<< "scale=6; ( ($stdev1E - $stdev2E) * 0.043 )" | xargs printf "%4.6f")
	echo "AvgStdDE $AvgStdDE $5" >> $dirpath/$1/AvgStdDE_$1.txt

	AvgDEprot=$(bc -l <<< "scale=6; ( ($pto1E - $pto2E) * 0.043 )" | xargs printf "%6.6f")
	echo "AvgDEprot $AvgDEprot $5">> $dirpath/$1/AvgDEprot_$1.txt
	AvgStdDEprot=$(bc -l <<< "scale=6; ( ($ptdev1E - $ptdev2E) * 0.043 )" | xargs printf "%6.6f")
	echo "AvgStdDEprot $AvgStdDEprot $5">> $dirpath/$1/AvgStdDEprot_$1.txt

}

while [[ $spMDjobidrp1status == RUNNING ]];
do
	sleep 5m
done


calverE $oxddir $reddir fe3crdfe3prm fe3crdfe2prm prm
calverE $reddir $oxddir fe2crdfe2prm fe2crdfe3prm prm

calverE $oxddir $reddir fe3crdfe3prm-rp fe3crdfe2prm-rp -rp
calverE $reddir $oxddir fe2crdfe2prm-rp fe2crdfe3prm-rp -rp

calverE $oxddir $reddir fe3crdfe3prm-rp1 fe3crdfe2prm-rp1 -rp1
calverE $reddir $oxddir fe2crdfe2prm-rp1 fe2crdfe3prm-rp1 -rp1

echo "Printing Marcus ET parameters for $protein $ligand complex " > $dirpath/MarcusET-parm.txt
echo "MD-traj AvgDEa AvgDEb lambdaET DGET lambdaETprot" >> $dirpath/MarcusET-parm.txt
calMarcusETparm () {
	echo "printing-for-trajectory $1"
	AvgDEred=$(cat $dirpath/$reddir/AvgDE_$reddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgDEoxd=$(cat $dirpath/$oxddir/AvgDE_$oxddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	lambdaET=$(bc -l <<< "scale=6; ((- $AvgDEred - $AvgDEoxd)/2)")
	DGET=$(bc -l <<< "scale=6; ((- $AvgDEred + $AvgDEoxd)/2)")

	AvgDEprotred=$(cat $dirpath/$reddir/AvgDEprot_$reddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	AvgDEprotoxd=$(cat $dirpath/$oxddir/AvgDEprot_$oxddir.txt | grep -- $1 | head -1 | awk '{print $2}')
	lambdaETprot=$(bc -l <<< "scale=6; ((- $AvgDEprotred - $AvgDEprotoxd)/2)")
	AvgDEred=$(bc -l <<< "scale=6; (- $AvgDEred)")
#	DGETprot=$(bc -l <<< "scale=6; (($AvgDEprotred + $AvgDEprotoxd)/2)"
	filename=$1
	if [[ $filename == prm ]]; 
	then
	        oxdcofacdis=$(tail -n +2 $dirpath/$oxddir/$protein-FMN-HEMcontacts.dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        redcofacdis=$(tail -n +2 $dirpath/$reddir/$protein-FMN-HEMcontacts.dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        avgcofacdis=$(bc -l <<< "scale=6; ( ( $oxdcofacdis + $redcofacdis )/2)" | xargs printf "%6.6f")
		echo "$1 $AvgDEoxd $AvgDEred $lambdaET $DGET $lambdaETprot $avgcofacdis" >> $dirpath/MarcusET-parm.txt
	else
		oxdcofacdis=$(tail -n +2 $dirpath/$oxddir/$protein-FMN-HEMcontacts"$1".dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        redcofacdis=$(tail -n +2 $dirpath/$reddir/$protein-FMN-HEMcontacts"$1".dat | awk '{print $4}' | st | tail -1 | awk '{print $5}')
	        avgcofacdis=$(bc -l <<< "scale=6; ( ( $oxdcofacdis + $redcofacdis )/2)" | xargs printf "%6.6f")
		echo "$1 $AvgDEoxd $AvgDEred $lambdaET $DGET $lambdaETprot $avgcofacdis" >> $dirpath/MarcusET-parm.txt
	fi
	
}

calMarcusETparm prm
calMarcusETparm -rp
calMarcusETparm -rp1
more $dirpath/$protein-$uniquename-MarcusET-parm.txt


