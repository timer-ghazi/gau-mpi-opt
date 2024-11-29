#!/bin/bash


#-------------------------------------------
DIM='\033[2m'   #Dim
BLD='\033[1m'   #Bold\brighter                                         |
DEF='\033[0m'   #Default color and effects
UND='\033[4m'   #Underline
CUR='\033[3m'   #Italic 
BIRed='\033[1;91m'        # Red

echo ""
echo    "                           -------------------------------"
echo -e "  ---" $BLD"Embarrassingly Parallel Gaussian Geometry Optimization with an External Code"$DEF "---" 
echo    "                      High-Level Energies + Gaussian Optimizer "
echo -e $DIM"                        Q.T., Marquette U., 2018--2022"$DEF
echo 
echo -e $BLD"                               Job Copy Utility"$DEF
echo 

if [ "$1"'t' == 't' ] ; then
echo -e " * How to use:"

echo -e     "   $BLD$(basename $0)$CUR old_file$DEF$BLD.com$CUR new_file"$DEF
echo -e $DIM"      - No extension for$CUR new_file$DEF"
echo -e $DIM"      - Any extension for$CUR old_file$DEF$DIM would work "
echo ""
exit
fi


if [ "$2"'t' == 't' ] ; then
echo -e " Please provide new name "
exit
fi

OLDN=$(echo $1 | sed 's/\.[^\.]*$//') 

if [ ! -f $OLDN.com ]
then
  echo -e " * Old input file $BLD$OLDN.com$DEF not found"
  exit
fi

if [ ! -f $OLDN.sh ]
then
  echo -e " * Old script file $BLD$OLDN.sh$DEF not found"
  exit
fi

NEWN=$2

OLDC="%OldChk="$OLDN.chk
NEWC="%Chk="$NEWN.chk

sed "s/%OldChk=.*$/$OLDC/;s/%Chk=.*$/$NEWC/" $OLDN.com > $NEWN.com 
echo ""
echo -e " * Copied information from the old input file $BLD$OLDN.com$DEF to the new one:  $BLD$NEWN.com$DEF "
echo -e "   - $OLDC "
echo -e "   - $NEWC "
echo ""
sed "s/$OLDN/$NEWN/" $OLDN.sh > $NEWN.sh
echo -e " * New Slurm submission file created: $BLD$NEWN.sh$DEF "
echo ""

echo -e $BLD" Next:"$DEF
echo -e "   - Edit the embedded Molpro/MRCC input in $BLD$NEWN.com"$DEF
echo -e "   - Submit the Slurm submission script: $BLD$NEWN.sh"$DEF
echo ""

if [ ! -f $OLDN.chk ]
then
  echo -e " $BIRed* Warning *$DEF Old checkpoint file $BLD$OLDN.chk$DEF not found"
  echo ""
fi

