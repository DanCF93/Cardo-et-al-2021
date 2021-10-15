#!/bin/bash

scriptName="1a_fastqc_pe_first"
jobName="1a"
jobMemory="4G"
jobLength="02:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}

sampleNames=($(${myDir}/read_param.sh sampleNames))
sampleNumber=${#sampleNames[@]}
fastqExtension=$(${myDir}/read_param.sh fastqExtension)
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaBatchQueue)
moduleFastqc=$(${myDir}/read_param.sh moduleFastqc)
moduleJava=$(${myDir}/read_param.sh moduleJava)

[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*
rm -f RUN/${scriptName}.sh 

for sample in "${sampleNames[@]}"
do
	[[ -d ${myDir}/../output/${sample}/fastqc/untrimmed/ ]] && rm -fr ${myDir}/../output/${sample}/fastqc/untrimmed/ ; mkdir -p ${myDir}/../output/${sample}/fastqc/untrimmed/
done

echo \#!/bin/bash > RUN/${scriptName}.sh
echo \#PBS -P ${project} >> RUN/${scriptName}.sh
echo \#PBS -q ${queue} >> RUN/${scriptName}.sh
echo \#PBS -N ${jobName} >> RUN/${scriptName}.sh
echo \#PBS -l select=1:ncpus=1:mem=${jobMemory} >> RUN/${scriptName}.sh
echo \#PBS -l walltime=${jobLength} >> RUN/${scriptName}.sh
echo \#PBS -o OUT/${scriptName}/ >> RUN/${scriptName}.sh
echo \#PBS -e ERR/${scriptName}/ >> RUN/${scriptName}.sh
echo \#PBS -J 1-${sampleNumber} >> RUN/${scriptName}.sh

echo module load ${moduleJava} >> RUN/${scriptName}.sh
echo module load ${moduleFastqc} >> RUN/${scriptName}.sh
echo sampleNames=\(${sampleNames[*]}\) >> RUN/${scriptName}.sh
echo fastqc -o ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/fastqc/untrimmed/ ${myDir}/../input/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\_1.${fastqExtension} ${myDir}/../input/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\_2.${fastqExtension} >> RUN/${scriptName}.sh
chmod u+x RUN/${scriptName}.sh 
qsub RUN/${scriptName}.sh

