#!/bin/bash

scriptName="4a_featurecount"
jobName="4a"
jobMemory="8G"
jobLength="05:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}

sampleNames=($(${myDir}/read_param.sh sampleNames))
sampleNumber=${#sampleNames[@]}
refGTF=$(${myDir}/read_param.sh refGTF)
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaBatchQueue)
moduleSAMTools=$(${myDir}/read_param.sh moduleSAMTools)
moduleFeatureCounts=$(${myDir}/read_param.sh moduleFeatureCounts)

[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*
rm -f RUN/${scriptName}.sh 

for sample in "${sampleNames[@]}"
do
	[[ -d ${myDir}/../output/${sample}/featurecount/ ]] && rm -fr ${myDir}/../output/${sample}/featurecount/ ; mkdir -p ${myDirr}../output/${sample}/featurecount/
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
echo module load ${moduleFeatureCounts} >> RUN/${scriptName}.sh

echo sampleNames=\(${sampleNames[*]}\) >> RUN/${scriptName}.sh

echo samtools sort -n ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.bam ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.sorted.markdup >> RUN/${scriptName}.sh
echo cd ${myDir}/../tmp \&\& featureCounts -O -p -F GTF -t exon -g gene_id -a ${myDir}/../resources/${refGTF} -o ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/featurecount/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.markdup.featurecount ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.sorted.markdup.bam >> RUN/${scriptName}.sh

echo samtools sort -n ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.bam ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.sorted.rmdup >> RUN/${scriptName}.sh
echo cd ${myDir}/../tmp \&\& featureCounts -O -p -F GTF -t exon -g gene_id -a ${myDir}/../resources/${refGTF} -o ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/featurecount/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.rmdup.featurecount ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.sorted.rmdup.bam >> RUN/${scriptName}.sh

chmod u+x RUN/${scriptName}.sh 
qsub RUN/${scriptName}.sh


