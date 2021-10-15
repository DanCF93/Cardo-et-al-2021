#!/bin/bash

scriptName="1b_index_genome"
jobName="1b"
jobMemory="40G"
jobCores="4"
jobLength="10:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}
refGenome=($(${myDir}/read_param.sh refGenome))
refGTF=($(${myDir}/read_param.sh refGTF))
seqLength=($(${myDir}/read_param.sh sequenceLength))
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaWorkQueue)
moduleSTAR=$(${myDir}/read_param.sh moduleSTAR)

[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*
rm -f RUN/${scriptName}.sh 

seqLength=`expr $seqLength - 1`

echo module load ${moduleSTAR} > RUN/${scriptName}.sh
echo STAR --runThreadN ${jobCores} --runMode genomeGenerate --genomeDir ${myDir}/../resources/ --genomeFastaFiles ${myDir}/../resources/$refGenome --sjdbGTFfile ${myDir}/../resources/${refGTF} --sjdbOverhang ${seqLength} >> RUN/${scriptName}.sh
chmod u+x RUN/${scriptName}.sh
qsub -P ${project} -q ${queue} -N ${jobName} -o OUT/${scriptName}.out -e ERR/${scriptName}.err -l select=1:ncpus=${jobCores}:mem=${jobMemory} -l walltime=${jobLength} RUN/${scriptName}.sh

