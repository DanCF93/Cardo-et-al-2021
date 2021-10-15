#!/bin/bash

scriptName="3a_mapping_qc"
jobName="3a"
jobMemory="12G"
jobLength="08:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}

sampleNames=($(${myDir}/read_param.sh sampleNames))
sampleNumber=${#sampleNames[@]}
refGenome=$(${myDir}/read_param.sh refGenome)
refGTF=$(${myDir}/read_param.sh refGTF)
paramIntersect=$(${myDir}/read_param.sh paramIntersect)
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaBatchQueue)
moduleJava=$(${myDir}/read_param.sh moduleJava)
modulePicard=$(${myDir}/read_param.sh modulePicard)
moduleSAMTools=$(${myDir}/read_param.sh moduleSAMTools)
moduleBEDTools=$(${myDir}/read_param.sh moduleBEDTools)
moduleBAMTools=$(${myDir}/read_param.sh moduleBAMTools)

[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*
rm -f RUN/${scriptName}.sh 

for sample in "${sampleNames[@]}"
do
	[[ -d ${myDir}/../output/${sample}/bamtools/ ]] && rm -fr ${myDir}/../output/${sample}/bamtools/ ; mkdir -p ${myDirr}../output/${sample}/bamtools/
	[[ -d ${myDir}/../output/${sample}/bedtools/ ]] && rm -fr ${myDir}/../output/${sample}/bedtools/ ; mkdir -p ${myDirr}../output/${sample}/bedtools/
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

echo module load ${modulePicard} >> RUN/${scriptName}.sh
echo module load ${moduleJava} >> RUN/${scriptName}.sh
echo module load ${moduleBEDTools} >> RUN/${scriptName}.sh
echo module load ${moduleSAMTools} >> RUN/${scriptName}.sh
echo module load ${moduleBAMTools} >> RUN/${scriptName}.sh

echo sampleNames=\(${sampleNames[*]}\) >> RUN/${scriptName}.sh

## mark duplicates and run bamtools stats

echo MarkDuplicates I=${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.multimap.Aligned.sortedByCoord.out.bam O=${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.bam M=${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.metrics.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT >> RUN/${scriptName}.sh
echo bamtools stats -in ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.bam \> ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/bamtools/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.stat.txt >> RUN/${scriptName}.sh
echo samtools index ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.bam >> RUN/${scriptName}.sh

## remove duplicates and run bamtools stats

echo MarkDuplicates I=${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.multimap.Aligned.sortedByCoord.out.bam O=${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.bam M=${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.metrics.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT >> RUN/${scriptName}.sh
echo bamtools stats -in ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.bam \> ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}//bamtools/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.stat.txt >> RUN/${scriptName}.sh
echo samtools index ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.bam >> RUN/${scriptName}.sh


chmod u+x RUN/${scriptName}.sh 
qsub RUN/${scriptName}.sh


