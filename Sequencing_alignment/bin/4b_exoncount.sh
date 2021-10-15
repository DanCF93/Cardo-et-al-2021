#!/bin/bash

scriptName="4b_exoncount"
jobName="4b"
jobMemory="8G"
jobLength="48:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}

sampleNames=($(${myDir}/read_param.sh sampleNames))
sampleNumber=${#sampleNames[@]}
refGTF=$(${myDir}/read_param.sh refGTF)
refGFF=${refGTF//.gtf/.gff}
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaBatchQueue)
moduleSAMTools=$(${myDir}/read_param.sh moduleSAMTools)
moduleHTSeq=$(${myDir}/read_param.sh moduleHTSeq)
moduleHTSeqPython=$(${myDir}/read_param.sh moduleHTSeqPython)


[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*
rm -f RUN/${scriptName}.sh 

for sample in "${sampleNames[@]}"
do
	[[ -d ${myDir}/../output/${sample}/exoncount/ ]] && rm -fr ${myDir}/../output/${sample}/exoncount/ ; mkdir -p ${myDirr}../output/${sample}/exoncount/
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

echo module load ${moduleSAMTools} >> RUN/${scriptName}.sh
echo module load ${moduleHTSeq} >> RUN/${scriptName}.sh
echo module load ${moduleHTSeqPython} >> RUN/${scriptName}.sh

echo sampleNames=\(${sampleNames[*]}\) >> RUN/${scriptName}.sh

echo samtools sort -n ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.bam ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.sorted >> RUN/${scriptName}.sh
echo samtools view ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.sorted.bam \| dexseq_count.py -p yes -a 10 -s no ${myDir}/../tmp/${refGFF} - ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/exoncount/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.exoncount.markdup >> RUN/${scriptName}.sh
#echo samtools view ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.sorted.bam \| sed \'s\/\$\/\\tNH:i:1\/g\' \| dexseq_count.py -p yes -a 10 -s no ${myDir}/../tmp/${refGFF} - ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/exoncount/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.exoncount.markdup >> RUN/${scriptName}.sh

echo samtools sort -n ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.bam ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.sorted >> RUN/${scriptName}.sh
echo samtools view ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.sorted.bam \| dexseq_count.py -p yes -a 10 -s no ${myDir}/../tmp/${refGFF} - ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/exoncount/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.exoncount.rmdup >> RUN/${scriptName}.sh
#echo samtools view ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.sorted.bam \| sed \'s\/\$\/\\tNH:i:1\/g\' \| dexseq_count.py -p yes -a 10 -s no ${myDir}/../tmp/${refGFF} - ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/exoncount/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.exoncount.rmdup >> RUN/${scriptName}.sh

chmod u+x RUN/${scriptName}.sh 
qsub RUN/${scriptName}.sh


