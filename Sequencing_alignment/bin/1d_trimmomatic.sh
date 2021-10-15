#!/bin/bash

scriptName="1d_trimmomatic"
jobName="1d"
jobMemory="32G"
jobLength="24:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}

sampleNames=($(${myDir}/read_param.sh sampleNames))
sampleNumber=${#sampleNames[@]}
fastqExtension=$(${myDir}/read_param.sh fastqExtension)
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaBatchQueue)
moduleTrimmomatic=$(${myDir}/read_param.sh moduleTrimmomatic)
moduleJava=$(${myDir}/read_param.sh moduleJava)

[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*
rm -f RUN/${scriptName}.sh 

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
echo module load ${moduleTrimmomatic} >> RUN/${scriptName}.sh
echo sampleNames=\(${sampleNames[*]}\) >> RUN/${scriptName}.sh

echo java -jar /software/genomics/trimmomatic/0.35/trimmomatic-0.35.jar PE -phred33 -trimlog ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed.log ${myDir}/../input/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\_1.${fastqExtension} ${myDir}/../input/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\_2.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_1.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.unpaired\_1.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_2.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.unpaired\_2.${fastqExtension} ILLUMINACLIP:${myDir}/../resources/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> RUN/${scriptName}.sh

chmod u+x RUN/${scriptName}.sh 
qsub RUN/${scriptName}.sh

