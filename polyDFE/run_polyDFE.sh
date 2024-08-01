#!/bin/bash
#PBS -N run_polyDFE
#PBS -l walltime=100:00:00
#PBS -l select=1:ncpus=1:mem=60gb:scratch_local=50gb
#PBS -j oe

if [ -z "$data" ]; then
        echo "Error! Specify relative path (relative to working directory) to the input file "
        echo "qsub -v 'data=polyDFEinput.txt'"
        exit 1
        fi



trap 'clean_scratch' TERM EXIT

if [ ! -d "$SCRATCHDIR" ] ; then echo "Scratch not created!" 1>&2; exit 1; fi 
bin="/storage/brno3-cerit/home/filip_kolar/programs/polyDFE/polyDFE"


cp $PBS_O_WORKDIR/$data $SCRATCHDIR/
cp /storage/brno3-cerit/home/filip_kolar/programs/polyDFE/input/init_A.txt $SCRATCHDIR/ # copy parameter files
cp /storage/brno3-cerit/home/filip_kolar/programs/polyDFE/input/init_C.txt $SCRATCHDIR/ # copy parameter files
cp /storage/brno3-cerit/home/filip_kolar/programs/polyDFE/input/params.basinhop $SCRATCHDIR/ # copy parameter files

data=`basename $data`

cd $SCRATCHDIR
### Full model A
###################
### 4 different models with different parameters estimated
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 2 -e > ${data/polydfe.input.txt/}_A_nor_noeps
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 3 -e > ${data/polydfe.input.txt/}_A_nor
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 4 -e > ${data/polydfe.input.txt/}_A_noeps
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 1 -e > ${data/polydfe.input.txt/}_A_full
### merge the results of the models
cat ${data/polydfe.input.txt/}_A_nor_noeps > ${data/polydfe.input.txt/}_A.all.models.txt
cat ${data/polydfe.input.txt/}_A_nor >> ${data/polydfe.input.txt/}_A.all.models.txt
cat ${data/polydfe.input.txt/}_A_noeps >> ${data/polydfe.input.txt/}_A.all.models.txt
cat ${data/polydfe.input.txt/}_A_full >> ${data/polydfe.input.txt/}_A.all.models.txt
###################
echo "model A finished at `date`"

### Full model C
##################
### again 4 different models
$bin -d $data -b params.basinhop 1 -m C -i init_C.txt 2 -e > ${data/polydfe.input.txt/}_C_nor_noeps
$bin -d $data -b params.basinhop 1 -m C -i init_C.txt 3 -e > ${data/polydfe.input.txt/}_C_nor
$bin -d $data -b params.basinhop 1 -m C -i init_C.txt 4 -e > ${data/polydfe.input.txt/}_C_noeps
$bin -d $data -b params.basinhop 1 -m C -i init_C.txt 1 -e > ${data/polydfe.input.txt/}_C_full # default

cat ${data/polydfe.input.txt/}_C_nor_noeps > ${data/polydfe.input.txt/}_C.all.models.txt
cat ${data/polydfe.input.txt/}_C_nor >> ${data/polydfe.input.txt/}_C.all.models.txt
cat ${data/polydfe.input.txt/}_C_noeps >> ${data/polydfe.input.txt/}_C.all.models.txt
cat ${data/polydfe.input.txt/}_C_full >> ${data/polydfe.input.txt/}_C.all.models.txt

echo "model C finished at `date`"


### Deleterious model

$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 6 -e > ${data/polydfe.input.txt/}_DEL_nor_noeps
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 7 -e > ${data/polydfe.input.txt/}_DEL_nor
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 8 -e > ${data/polydfe.input.txt/}_DEL_noeps
$bin -d $data -b params.basinhop 1 -m A -i init_A.txt 5 -e > ${data/polydfe.input.txt/}_DEL_full

cat ${data/polydfe.input.txt/}_DEL_nor_noeps > ${data/polydfe.input.txt/}_DEL.all.models.txt
cat ${data/polydfe.input.txt/}_DEL_nor >> ${data/polydfe.input.txt/}_DEL.all.models.txt
cat ${data/polydfe.input.txt/}_DEL_noeps >> ${data/polydfe.input.txt/}_DEL.all.models.txt
cat ${data/polydfe.input.txt/}_DEL_full >> ${data/polydfe.input.txt/}_DEL.all.models.txt

echo "model DEL finished at `date`"


cp $SCRATCHDIR/* ${PBS_O_WORKDIR}/ || export CLEAN_SCRATCH=false

