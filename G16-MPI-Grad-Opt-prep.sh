#!/bin/bash

CHARGE=0
MULT=1

ESTCODE=MOLPRO
#ESTCODE=MRCC

# how many cores will each MRCC run use
MRCCCORES=14
MRCCMEM=500GB

# how many cores will each Molpro run use
MOLPROCORES=6
# memory per core in gigawords 
MOLPROMEM=3

# Scratch dir on the nodes 
LOCAL_SCRATCH_DIR="/local"

#############################################################
#
#  Scratch disk free space check 
# 
#  If enabled, will first check how much free space is available 
#  on the assigned nodes, and if there is less than 
#  $MIN_SCRATCH_SPACE_PER_SINGLEPOINT GB per single-point, 
#  will stop the job.
#  Also prints SLURM exclude directive to avoid nodes with insufficent space 
# 
CHECK_SCRATCH_SPACE="YES"
MIN_SCRATCH_SPACE_PER_SINGLEPOINT=550
#
#############################################################

PYSCRDIR=/mmfs1/home/timerghazinq/local/gau-opt
CONDADIR=/mmfs1/home/timerghazinq/anaconda3
MRCCDIR=/mmfs1/home/timerghazinq/local/MRCC2022
MOLPRODIR=/mmfs1/home/timerghazinq/local/molpro/molprop_2015_1_linux_x86_64_i8/bin

cat > .mrcci << EOF
basis=cc-pVDZ
calc=CCSD
mult=$MULT
charge=$CHARGE
mem=$MRCCMEM
geom=xyz
%s
EOF

cat > .molpro << EOF
***,SinglePt
memory,$MOLPROMEM,G
gthresh,energy=1.d-8,orbital=1.d-8,grid=1.d-8
!NOSYM;
NOORIENT;
gparam,dsyevd=0
geomtyp=xyz
! do not modify--XYZ geometry will be inserted 
GEOMETRY={%s
}
BASIS
default=vdz
END
set,charge=$CHARGE
hf;               
show,Energy
!  the names are set by the script; do not modify:
put,xyz,%s.xyz;
EOF


#-------------------------------------------
DIM='\033[2m'   #Dim
BLD='\033[1m'   #Bold\brighter                                         |
DEF='\033[0m'   #Default color and effects
UND='\e[4m'   #Underline
CUR='\e[3m'   #Italic font

echo ""
echo    "                           -------------------------------"
echo -e "  ---" $BLD"Embarrassingly Parallel Gaussian Geometry Optimization with an External Code"$DEF "---" 
echo    "                      High-Level Energies + Gaussian Optimizer "
echo    "                          Active Interface: $ESTCODE "
echo -e $DIM"                        Q.T., Marquette U., 2018--2024"$DEF
echo -e $DIM"                           qadir.timerghazin@mu.edu"$DEF
echo ""

if [ "$1"'t' == 't' ] ; then

#echo -e $BLD" * Currently supported electronic structure theory (EST) codes:"$DEF
echo -e " * Currently supported electronic structure theory (EST) codes:"$DEF
echo -e "    - "$BLD"MRCC"$DEF
echo -e "    - "$BLD"Molpro"$DEF
echo ""


echo -e $BLD" * How to use:"$DEF

echo -e "   0. Make sure to copy this script into your work directory:"
echo -e "    $BLD  cp $PYSCRDIR/$(basename $0) . "$DEF
echo -e $DIM"      The script may require minimal customization such as molecular charge and"$DEF
echo -e $DIM"      multiplicity on top of the file."$DEF

echo -e "   1. Run this script with a XYZ file name: $BLD./$(basename $0)$CUR mymolecule$DEF$BLD.xyz"$DEF
echo -e $DIM"      This will rotate the coordinates into standard orientation and also"
echo -e     "      generate$CUR mymolecule$DEF$DIM.symm file with proposed symmetry constraints"$DEF

echo -e "   2. Correct any issues with generated constraints in$BLD$CUR mymolecule$DEF$BLD.symm"$DEF
echo -e $DIM"      E.g. any coordinate equal to 0.0 is assumed to be in a symmetry plane" 
echo -e     "      with zero force ('X00 = ZERO'). The algorithm is primitive--human editing is required!"
echo -e     "      Make sure to remove any constraints for "$CUR"C"$DEF$DIM"₁ molecules!"$DEF

echo -e "   3. Run this script again, now with two files: "
echo -e "      $BLD./$(basename $0)$CUR mymolecule$DEF$BLD.xyz$CUR mymolecule$DEF$BLD.symm$CUR$DIM optional_job_name"$DEF 
echo -e $DIM"      This will generate necessary input files and scripts setup for your system."$DEF
echo -e $DIM"      If specified,$CUR optional_job_name$DEF$DIM will be used as the main file name."$DEF
echo -e $DIM"       - To daisy-chain several optimizations (e.g. with increasing level of theory) use:"$DEF
echo -e $DIM"         $BLD./$(basename $0)$CUR mymolecule$DEF$DIM$BLD.xyz$CUR mymolecule$DEF$BLD$DIM.symm$CUR$DIM job1 job2 job3 ..."$DEF
echo -e $DIM"         This will create$CUR job1$DEF$DIM.com,$CUR job2$DEF$DIM.com, etc. files sequentially run from the "$DEF
echo -e $DIM"        $CUR mymolecule$DEF$DIM.sh script, each job starting with the previous job's checkpoint file."$DEF

echo -e "   4. Run preliminary DFT frequency job$BLD dft-"$CUR"mymolecule$DEF$BLD-Freq.com$DEF as a regular Gaussian job"
echo -e $DIM"      This will provide a good starting Hessian for the geometry optimization."$DEF

echo -e "   5. Edit "$BLD$CUR"mymolecule$DEF$BLD.com$DEF to setup the $ESTCODE calculation and Gaussian optimization"
echo -e $DIM"      The $ESTCODE input header is embedded in the Gaussian .com file as a comment, so "$DEF
echo -e $DIM"      each line starts with '!' that will be automatically removed later. "$DEF
echo -e $DIM"      Also, the $ESTCODE input header is defined on top of this script (no '!' there). "$DEF

echo -e "   6. Adjust, if needed, and submit the Slurm job script "$BLD$CUR"mymolecule$DEF$BLD.sh$DEF"
echo -e $DIM"      - if at all possible, run only one MRCC job per node ('#SBATCH --ntasks-per-node=1') "$DEF
echo -e $DIM"      - monitor the messages in a Slurm output file$CUR mymolecule$DEF$DIM-*.SLURM for job progress/errors "$DEF
echo ""

echo -e " " $DIM$UND"Note:"$DEF
echo -e $DIM"  * Only simple symmetry constraints implemented, but they can significantly cut the number of          "$DEF
echo -e $DIM"    individual $ESTCODE energy evaluations; automatically generated constraints "$UND"will$DEF$DIM have mistakes,  "$DEF
echo -e $DIM"    so editing the .symm file is always necessary. "$DEF
echo -e $DIM"    For "$CUR"C"$DEF$DIM"₁ molecules, make sure that there are no constraints (i.e. 'X01 = X01', etc. for all)          "$DEF
echo ""


echo -e $DIM"  * It is recommended to start with a lower coupled-cluster level optimization, and use the optimized     "$DEF 
echo -e $DIM"    geometry (and force constants) as a starting point for optimization with the next, higher level,      "$DEF 
echo -e $DIM"    which in turn will be the starting point for the higher yet theory-level optimization.                "$DEF 
echo -e $DIM"    E.g., to perform CCSDT(Q) optimization, a CCSD-->CCSD(T)-->CCSDT(Q) sequence can be used.             "$DEF 
echo -e $DIM"    Two options to do this:                                                                               "$DEF 
echo -e $DIM"     - After an optimization is complete, use$BLD copy2next-Gau-MRCC.sh$DEF$DIM script to produce the "$DEF 
echo -e $DIM"       input & submittion script for the next step (.chk file from the previous job will be used)         "$DEF 
echo -e $DIM"     - Generate and run the entire sequence at once (see above):                                          "$DEF 
echo -e $DIM"       $BLD./$(basename $0)$CUR mymolecule$DEF$DIM$BLD.xyz$CUR mymolecule$DEF$BLD$DIM.symm$CUR mymolecule-CCSD mymolecule-CCSDt mymolecule-CCSDTq "$DEF 

echo ""

echo -e "  "$DIM$UND"Requires:"$DEF
echo -e $DIM"   * Rest of this package--Python scripts in $BLD\$PYSCRDIR"$DEF
echo -e $DIM"   * MRCC (in $BLD\$MRCCDIR"$DEF$DIM") and Gaussian (as a module)"$DEF 
echo -e $DIM"   * Anaconda in $BLD\$CONDADIR"$DEF$DIM "with the following packages installed:"$DEF
echo -e $DIM"      - Python 2.7 (installed as 'py2' environment)"$DEF
echo -e $DIM"      - numpy"$DEF
echo -e $DIM"      - mpi4py (which also installs mpich MPI library)"$DEF

echo ""
echo -e $BLD"  * Current settings (edit this script to correct):"$DEF
echo -e "    - Molecule charge: $BLD$CHARGE$DEF multiplicity: $BLD$MULT$DEF "
echo -e "    - Active EST code interface: $BLD$ESTCODE$DEF"
echo -e "    - MRCC OpenMP cores: $BLD$MRCCCORES$DEF, MRCC memory: $BLD$MRCCMEM$DEF"
echo -e "    - Molpro cores: $BLD$MOLPROCORES$DEF, Molpro memory/core, GW: $BLD$MOLPROMEM$DEF"

echo ""
echo -e $BLD"  * Current directory settings (edit this script to correct):"$DEF
echo -e "    - Python scripts: PYSCRDIR=$BLD$PYSCRDIR$DEF" 
echo -e "    - Anaconda: CONDADIR=$BLD$CONDADIR$DEF" 
echo -e "    - MRCC: MRCCDIR=$BLD$MRCCDIR$DEF" 
echo -e "    - Molpro: MOLPRODIR=$BLD$MOLPRODIR$DEF" 
echo ""

rm .mrcci
rm .molpro
exit

fi

source $CONDADIR/etc/profile.d/conda.sh
conda activate py2 
PATH=$PYSCRDIR:$PATH
PYTHONPATH=$PYSCRDIR:$PYTHONPATH


if [ "$2"'t' == 't' ] ; then
#-----------------------------------------------------------------------------------
# Producing Standard Orientaton 
#
XYZF=$1
if [ ! -f $XYZF ]
then
   echo $XYZF "not found"
   exit 
fi

echo -e $BLD" Job Setup Part 1:$DEF Rotating the XYZ file into standard orientation and "
echo -e "                   generating symmetry constraints template"
echo ""

echo -e $BLD" * Original XYZ file $XYZF:"$DEF
cat $XYZF

COORDS=$(awk  '{ if ((NR != 2) && (NR != 1)) print $0}' $XYZF )

COMF=tmp-${1/xyz/com}
LOGF=tmp-${1/xyz/log}

cat > $COMF << EOF
#P guess=only  
Symmetry(Loose,PrintOrientation)


$CHARGE $MULT
$COORDS


EOF

echo -e $DIM
echo "Executing g16 to produce standard orientation..."
module load gaussian > _err 2>&1
cat _err
rm _err
g16 $COMF 
echo -e $DEF

awk 'BEGIN{aS[-1]="X";aS[0]="X";aS[1]="H";aS[2]="He";aS[3]="Li";aS[4]="Be";aS[5]="B";aS[6]="C";
aS[7]="N";aS[8]="O";aS[9]="F";aS[10]="Ne";aS[11]="Na";aS[12]="Mg";aS[13]="Al";aS[14]="Si";aS[15]="P";
aS[16]="S";aS[17]="Cl";aS[18]="Ar";aS[19]="K";aS[20]="Ca";aS[21]="Sc";aS[22]="Ti";aS[23]="V";aS[24]="Cr";
aS[25]="Mn";aS[26]="Fe";aS[27]="Co";aS[28]="Ni";aS[29]="Cu";aS[30]="Zn";aS[31]="Ga";aS[32]="Ge";aS[33]="As";
aS[34]="Se";aS[35]="Br";aS[36]="Kr";aS[46]="Pd";aS[47]="Ag";aS[48]="Cd";aS[50]="Sn";aS[51]="Sb";aS[53]="I";
aS[54]="Xe";aS[64]="Pt";aS[65]="Au";aS[66]="Hg";aS[67]="Tl";aS[68]="Pb";aS[69]="Bi"}
{if ($1=="Coordinates" && $2=="in" && $3=="use:") {getline; i=0; while ($1 != "Leave") {i++; at[i]=$1; x[i]=$2; y[i]=$3; z[i]=$4; getline} } }
END{print i; print "Standard orientation";for (k=1;k<=i;k++){printf("%s   % 16.12f   % 16.12f   % 16.12f\n",aS[at[k]],x[k],y[k],z[k])}}' $LOGF > std_

mv $XYZF original.$XYZF
mv std_ $XYZF

echo -e $DIM"Original file saved as original.$XYZF"$DEF
echo ""
echo -e $BLD" * Standard orientation generated:"$DEF
cat $XYZF

echo ""

echo -e $BLD" * Generating symmetry constraint template/proposal:"$DEF
gradienteval.py $XYZF --batch
PROPSYM=${XYZF/.xyz/.symm}

echo -e " * Proposed symmetry constraints for editing saved in $BLD$PROPSYM$DEF file. "
echo ""

echo -e $BLD" Next:"$DEF
echo -e "   - Fix the generated symmetry constraints in $BLD$PROPSYM$DEF "
echo -e "   - Run this script again: $BLD./$(basename $0) $XYZF $PROPSYM $DEF"
echo ""

rm .mrcci

# 
#  PART 1 END
#
#-----------------------------------------------------------------------------------
else
#-----------------------------------------------------------------------------------
# 
#  PART 2 BEGIN
#
# if we have symmetry contstraints, we're good to go
XYZF=$1
SYMF=$2

if [ ! -f $XYZF ]
then
   echo $XYZF "not found"
   exit 
fi
if [ ! -f $SYMF ]
then
   echo $SYMF "not found"
   exit 
fi

# cleaning up after step #1:
if [ -f original.$XYZF ]
then
   rm original.$XYZF
fi
if [ -f tmp-${1/xyz/com} ]
then
    rm tmp-${1/xyz/com}
fi
if [ -f tmp-${1/xyz/log} ]
then
    rm tmp-${1/xyz/log}
fi

if [ $ESTCODE == 'MRCC' ] ; then
    ESTCORES=$MRCCCORES
    ESTMEM=$MRCCMEM
elif [ $ESTCODE == 'MOLPRO' ] ; then
    ESTCORES=$MOLPROCORES
    ESTMEM=$(( $MOLPROMEM*8*$MOLPROCORES ))GB
fi

echo -e $BLD" Job Setup Part 2: Parsing Symmetry Constraints and Deploying Files" $DEF
echo ""

BASEN=${1/.xyz/}

gradienteval.py $XYZF -s $SYMF
# number of instances 
NINSTANCES=$(gradienteval.py $XYZF -s $SYMF --batch)
DFTCOM=dft-$BASEN-Freq.com 

COORDS=$(awk  '{ if ((NR != 2) && (NR != 1)) print $0}' $XYZF )
SYMMETRYCONSTR=$(cat $SYMF)

echo ""
echo -e $BLD" Setting up calculation parameters:" $DEF
echo -e " * $BLD$NINSTANCES$DEF parallel $ESTCODE energy evaluations, each using $BLD$ESTCORES$DEF cores"


echo ""
echo -e $BLD" Deploying main job files (may require user modifications):" $DEF

cat > $DFTCOM << EOF
--Link1--
%Nproc=8
%Mem=8GB
%Chk=$BASEN-DFT-Freq.chk
#P wB97XD/Def2TZVPP
Freq=noraman

DFT run to generate a good initial Hessian for External Code+Gaussian optimization for $BASEN

$CHARGE $MULT
$COORDS

EOF
echo -e " * $BLD$DFTCOM$DEF: Gaussian input to generate initial Hessian with DFT (run it first!)"



if [ "$3"'t' == 't' ] ; then
   # only 2 parameters, single .com file with default name 
   STEPS=$BASEN
else
   # at least one step name 
   STEPS=${@:3}
fi   

# script header  
SLURMNAME=$BASEN
SSCR=$BASEN.sh
if [ "$4"'t' == 't' ] ; then
    if [ "$3"'t' != 't' ] ; then
        SLURMNAME=$3
        SSCR=$3.sh
    fi
fi

cat > $SSCR << EOF
#!/bin/sh
#SBATCH --job-name=$SLURMNAME
#SBATCH --output=%x-%j.SLURM
#SBATCH --ntasks=$NINSTANCES
#SBATCH --cpus-per-task=$ESTCORES
#SBATCH --mem=$ESTMEM
# #SBATCH --ntasks-per-node=1
# * Note: use --ntasks-per-node=1 if at all possible--MRCC runs seem to significantly interfere with each other 

#  #SBATCH --partition=ai
#  #SBATCH --nodes=1
#  #SBATCH --exclude=raj-gpu09
#  #SBATCH --exclude=raj-cn[12-24,26-29]
#  #SBATCH --dependency=afterok:11254323

INSTANCES=$NINSTANCES


export WKDIR=\$(pwd)

module load gaussian

source $CONDADIR/etc/profile.d/conda.sh
conda activate py2
export PATH=$PYSCRDIR:\$PATH
export PYTHONPATH=$PYSCRDIR:\$PYTHONPATH

MY_NODEFILE=/tmp/hosts.\$SLURM_JOB_ID
scontrol show hostnames > \$MY_NODEFILE

EOF

#############################################################
#
#  Scratch disk free space check 
# 
#  If enabled, will first check how much free space is available 
#  on the assigned nodes, and if there is less than 
#  $MIN_SCRATCH_SPACE_PER_SINGLEPOINT GB per single-point, 
#  will stop the job.
#  Also prints SLURM exclude directive to avoid nodes with insufficent space 
# 

cat >> $SSCR << EOF

REQUIRED_SPACE_GB=$MIN_SCRATCH_SPACE_PER_SINGLEPOINT

# Function to check available space on a node
check_space() {
    local node=\$1
    local tasks_per_node=\$2
    local required_space=\$((tasks_per_node * REQUIRED_SPACE_GB))
    
    local df_output=\$(ssh \$node "df -BG $LOCAL_SCRATCH_DIR | tail -1")
    total_gb=\$(echo "\$df_output" | awk '{print \$2}' | sed 's/G//')
    space_gb=\$(echo "\$df_output" | awk '{print \$4}' | sed 's/G//')
    space_percent=\$(echo "\$df_output" | awk '{print \$5}' | sed 's/%//')

    sufficient_marker="YES"
    if [ \$space_gb -lt \$required_space ]; then
         sufficient_marker=" NO"
    fi
    
    printf " %-20s %10s %10s %10s %10s %10s\n" "\$node" "\$total_gb GB" "\$space_gb GB" "\$space_percent%" "\$required_space GB" "\$sufficient_marker"
    
    if [ \$space_gb -lt \$required_space ]; then
        return 1
    fi
    
    return 0
}


# Check space on each node
num_nodes=\$(cat \$MY_NODEFILE | wc -l)
tasks_per_node=\$((INSTANCES / num_nodes))
remaining_tasks=\$((INSTANCES % num_nodes))
insufficient_nodes=()

echo "Nodes assigned:"
cat \$MY_NODEFILE
echo ""
echo "INSTANCES=\$INSTANCES"
echo "num_nodes=\$num_nodes"
echo ""

echo "---------------------------------------------------------------------------------------"
printf " % -20s %10s %10s %10s %10s %10s\n" "Node" "Capacity" "Free" "Free" "Required" "Enough?" 
printf " % -20s %10s %10s %10s %10s\n" "" "(GB)" "(GB)" "(%)" "(GB)" 
echo "---------------------------------------------------------------------------------------"
for node in \$(cat \$MY_NODEFILE); do
    if [ \$remaining_tasks -gt 0 ]; then
        check_space \$node \$((tasks_per_node + 1))
        check_status=\$?
        remaining_tasks=\$((remaining_tasks - 1))
    else
        check_space \$node \$tasks_per_node
        check_status=\$?
    fi
    
    if [ \$check_status -eq 1 ]; then
        insufficient_nodes+=(\$node)
    fi
done
echo "---------------------------------------------------------------------------------------"

if [ \${#insufficient_nodes[@]} -eq 0 ]; then
    echo "All nodes have sufficient space. All clear!"
else
    echo ""
    echo "The following nodes have insufficient space:"
    printf '%s\n' "\${insufficient_nodes[@]}"

    exclude_directive="#SBATCH --exclude=\$(IFS=,; echo "\${insufficient_nodes[*]}")"
    echo ""
    echo "To exclude the insufficient nodes, add the following directive to your script:"
    echo "\$exclude_directive"

    exit 1
fi

EOF
#
#
#############################################################

OLDCHK=$BASEN-DFT-Freq.chk
COMLST=""

if [ "$4"'t' != 't' ] ; then
    echo -e " * Gaussian input files with embedded $ESTCODE input and symmetry constraints:"
fi

# ***********************************************************************************************
for ASTEP in $STEPS
do

CHK=$ASTEP.chk
COM=$ASTEP.com 
COMLST=$COMLST" "$COM

if [ $ESTCODE == 'MRCC' ] ; then
    INPTEMPL=$( awk '{print "!"$0}' .mrcci )
elif [ $ESTCODE == 'MOLPRO' ] ; then
    INPTEMPL=$( awk '{print "!"$0}' .molpro )
fi

cat > $COM << EOF
!CODE $ESTCODE
!# input file template
!****** $ESTCODE TEMPLATE BEGIN
$INPTEMPL
!****** $ESTCODE TEMPLATE END
!****** SYMMETRY CONSTRAINTS BEGIN
$SYMMETRYCONSTR
!****** SYMMETRY CONSTRAINTS END
--Link1--
%Chk=$CHK
%OldChk=$OLDCHK
#P External='$PYSCRDIR/gau-external.py'
Opt(Tight,ReadFC,NoMicro,MaxCycle=15)
geom=check

$ESTCODE+Gaussian optimization for $BASEN

$CHARGE $MULT

EOF

rm .mrcci
rm .molpro

if [ "$4"'t' == 't' ] ; then
    echo -e " * $BLD$COM$DEF: Gaussian input file with embedded $ESTCODE input and symmetry constraints"
else
    #echo -e "   - $BLD$COM$DEF $DIM(%OldChk=$OLDCHK, %Chk=$CHK)$DEF"
    printf  "   - $BLD%-15s$DEF\t$DIM\t%%OldChk=%20s\t%%Chk=%20s$DEF\n" $COM $OLDCHK $CHK
fi


cat >> $SSCR << EOF

# if an EST run ends with an error, the output files are collected in ERRORS/
rm -rf ERRORS &> /dev/null

JOBFILE=$ASTEP
echo "Starting MPI for \$JOBFILE on \$(date) Host: \$(hostname)" >> \$JOBFILE.times
echo "Starting MPI for \$JOBFILE on \$(date) Host: \$(hostname)" 
echo "Hosts used: \$(cat \$MY_NODEFILE)"
mpirun -np \$INSTANCES -f \$MY_NODEFILE -launcher ssh  gau-mpigrad-opt.py \$JOBFILE.com 



cat \$JOBFILE-?0?-*.time >> \$JOBFILE.times
tar czf \$JOBFILE.tgz \$JOBFILE-?0?-* \$JOBFILE-MOLECULE.*

mv \$JOBFILE-MOLECULE.out    \$JOBFILE.out  &> /dev/null
mv \$JOBFILE-MOLECULE.MINP   \$JOBFILE.MINP &> /dev/null
mv \$JOBFILE-MOLECULE.mpi   \$JOBFILE.mpi &> /dev/null
mv \$JOBFILE-MOLECULE.inp   \$JOBFILE.inp &> /dev/null
rm \$JOBFILE-?0?-* \$JOBFILE-MOLECULE.*

rm gau_EOut &> /dev/null  
rm GAU_DONE &> /dev/null

EOF
OLDCHK=$CHK

done 
# ***********************************************************************************************

if [ "$4"'t' == 't' ] ; then
    echo -e " * $BLD$SSCR$DEF: Slurm submission file for the main $ESTCODE/Gaussian run"
else
    echo -e " * $BLD$SSCR$DEF: Slurm submission file for the main $ESTCODE/Gaussian run, .com files will be run sequentially" 
fi
echo ""

###############################################################################################
echo -e  $DIM $UND"Deploying auxillary scripts"$DEF

#----------------------------------------------------------------------------
#       __  __ ____   ____ ____   ____            _       _          
#      |  \/  |  _ \ / ___/ ___| / ___|  ___ _ __(_)_ __ | |_        
#      | |\/| | |_) | |  | |     \___ \ / __| '__| | '_ \| __|       
#      | |  | |  _ <| |__| |___   ___) | (__| |  | | |_) | |_        
#      |_|  |_|_| \_\\____\____| |____/ \___|_|  |_| .__/ \__|       
#                                                  |_|               
                                                                     
if [ $ESTCODE == 'MRCC' ] ; then

cat > run-mrcc.sh << EOF
#!/bin/bash

BASEN=\${1/.MINP/}

export SCRDIR=\$(mktemp -d -p $LOCAL_SCRATCH_DIR)
export WKDIR=\$(pwd)

cd \$SCRDIR

export OMP_NUM_THREADS=$MRCCCORES
export MKL_NUM_THREADS=$MRCCCORES

CPUVENDOR=\$(lscpu | grep "Vendor ID" | awk '{print \$3}')

if [ "\$CPUVENDOR" == "AuthenticAMD" ]; then
    export KMP_AFFINITY=compact
else 
    export KMP_AFFINITY=scatter
fi

HST=\$(hostname)
echo "\$(date) Running \$1 on \$HST Scratch: \$SCRDIR" > \$WKDIR/\$BASEN.messages
echo "KMP_AFFINITY="\$KMP_AFFINITY                 >> \$WKDIR/\$BASEN.messages
echo "OMP_NUM_THREADS="\$OMP_NUM_THREADS           >> \$WKDIR/\$BASEN.messages
echo "MKL_NUM_THREADS="\$MKL_NUM_THREADS           >> \$WKDIR/\$BASEN.messages


if [ -f \$WKDIR/GENBAS ]
then
   cp \$WKDIR/GENBAS .     
fi

cp \$WKDIR/\$1 MINP 

export PATH=$MRCCDIR:\$PATH 

/usr/bin/time -f "%e %P" -o \$BASEN.time dmrcc > \$WKDIR/\$BASEN.out


EXSTATUS=\$?

# Check the exit status
if [ \$EXSTATUS -ne 0 ]; then
  echo "\$BASEN The executable terminated with an error \$EXSTATUS" >> \$WKDIR/\$BASEN.messages
  mkdir -p \$WKDIR/ERRORS && cp \$WKDIR/\$BASEN.* \$WKDIR/ERRORS/
else
  echo "\$BASEN \$(cat \$BASEN.time) Finished on: \$(date) Host: \$HST" >> \$WKDIR/\$BASEN.time 
fi

rm -rf \$TMPDIR
cd \$WKDIR 
rm -rf \$SCRDIR

# Exit the script with the same status as the executable
exit \$EXSTATUS


EOF
chmod +x run-mrcc.sh
echo -e $DIM "*" $UND"run-mrcc.sh$DEF$DIM: MRCC execution script$DEF"

#---------------------------------------------------------------------------------
#  __  __       _                   ____            _       _
# |  \/  | ___ | |_ __  _ __ ___   / ___|  ___ _ __(_)_ __ | |_
# | |\/| |/ _ \| | '_ \| '__/ _ \  \___ \ / __| '__| | '_ \| __|
# | |  | | (_) | | |_) | | | (_) |  ___) | (__| |  | | |_) | |_
# |_|  |_|\___/|_| .__/|_|  \___/  |____/ \___|_|  |_| .__/ \__|
#                |_|                                 |_|

elif [ $ESTCODE == 'MOLPRO' ] ; then

cat > run-molpro.sh << EOF
#!/bin/bash

BASEN=\${1/.mpi/}

export MOLPRO_SCRDIR=\$(mktemp -d -p $LOCAL_SCRATCH_DIR)
export TMPDIR=\$MOLPRO_SCRDIR
WKDIR=\$(pwd)

MLPCPUS=$MOLPROCORES

HST=\$(hostname)
echo "\$(date) Running \$1 on \$HST Scratch: \$TMPDIR" > \$WKDIR/\$BASEN.messages
echo "Number of CPUs used "\$MLPCPUS >> \$WKDIR/\$BASEN.messages

module unload slurm 

hydra_exe=$MOLPRODIR/mpiexec.hydra
molpro_exe=$MOLPRODIR/molpro.exe

/usr/bin/time -f "%e %P" -o \$WKDIR/\$BASEN.time \$hydra_exe -launcher ssh -hosts \$HST -np \$MLPCPUS \$molpro_exe \$1 -t 1  &>> \$WKDIR/\$BASEN.messages

EXSTATUS=\$?

# Check the exit status
if [ \$EXSTATUS -ne 0 ]; then
  echo "\$BASEN The executable terminated with an error \$EXSTATUS" >> \$WKDIR/\$BASEN.messages
  mkdir -p \$WKDIR/ERRORS && cp \$WKDIR/\$BASEN.* \$WKDIR/ERRORS/
fi

echo "\$BASEN \$(cat \$BASEN.time) Finished on: \$(date) Host: \$HST" >> \$WKDIR/\$BASEN.time 

rm -rf \$TMPDIR

# Exit the script with the same status as the executable
exit \$EXSTATUS

EOF
chmod +x run-molpro.sh
echo -e $DIM "*" $UND"run-molpro.sh$DEF$DIM: Molpro execution script$DEF"

fi

#######################################################################################################################3



cat > gaussian.sh << EOF
#!/bin/bash
source $CONDADIR/etc/profile.d/conda.sh
conda activate py2
export PATH=$PYSCRDIR:\$PATH
export PYTHONPATH=$PYSCRDIR:\$PYTHONPATH

module load gaussian
export GAUSS_SCRDIR=$PWD

echo "gaussian.sh started at \$(hostname)"
echo "Parameter taken by gaussian.sh:" \$1

g16 \$1

EXSTATUS=\$?

# Exit the script with the same status as the executable
exit \$EXSTATUS

EOF
chmod +x gaussian.sh 
echo -e $DIM "*" $UND"gaussian.sh$DEF$DIM: Gaussian wrapper script$DEF"
echo ""

echo -e $BLD" Next:"$DEF
echo -e "   - Run a preliminary DFT frequency calculation:$BLD qg $DFTCOM"$DEF 

if [ "$4"'t' == 't' ] ; then
    echo -e "   - Edit embedded $ESTCODE input in$BLD$COMLST$DEF"
    echo -e "   - Adjust $BLD$SSCR$DEF script as needed" 
    echo -e "   - Submit the script:$BLD sbatch $SSCR$DEF" 
    echo ""
    echo -e "   - To continue optimization using a higher level of theory, use$BLD copy2next-Gau-MRCC.sh$DEF script " 
    echo -e $DIM"     • It is recommended to start with a lower coupled-cluster level optimization, e.g. CCSD(T),           "$DEF
    echo -e $DIM"       then use$BLD copy2next-Gau-MRCC.sh$DEF$DIM to produce input & submittion script for the next step,  "$DEF
    echo -e $DIM"       e.g. CCSDT(Q); the new job thus generated will start off with the geometry and                      "$DEF
    echo -e $DIM"       force constants of the first one.                                                                   "$DEF
else
    echo -e "   - Edit embedded MRCC inputs in each .com file:$BLD$COMLST$DEF"
    echo -e "   - Adjust $BLD$SSCR$DEF script as needed" 
    echo -e "   - Submit the script:$BLD sbatch $SSCR$DEF"
fi



echo ""
#-------------------------------------------------------------------
fi






