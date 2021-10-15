#!/bin/bash

scriptName="2b_star_align"
jobName="2b"
jobMemory="40G"
jobCores="1"
jobLength="24:00:00"

myDir=$(pwd)
tmp=$( echo $myDir | sed "s/.*an0\([0-9]\+\).*/\1/p" -n )
jobName=${jobName}.${tmp}

sampleNames=($(${myDir}/read_param.sh sampleNames))
sampleNumber=${#sampleNames[@]}
refGenome=$(${myDir}/read_param.sh refGenome)
fastqExtension=$(${myDir}/read_param.sh fastqExtension)
project=$(${myDir}/read_param.sh arccaProject)
queue=$(${myDir}/read_param.sh arccaBatchQueue)
moduleSTAR=$(${myDir}/read_param.sh moduleSTAR)

[ -d OUT/${scriptName}/ ] || mkdir OUT/${scriptName}/ ; rm -f OUT/${scriptName}/*
[ -d ERR/${scriptName}/ ] || mkdir ERR/${scriptName}/ ; rm -f ERR/${scriptName}/*

rm -f RUN/${scriptName}.sh 

echo \#!/bin/bash > RUN/${scriptName}.sh
echo \#PBS -P ${project} >> RUN/${scriptName}.sh
echo \#PBS -q ${queue} >> RUN/${scriptName}.sh
echo \#PBS -N ${jobName} >> RUN/${scriptName}.sh
echo \#PBS -l select=${jobCores}:ncpus=${jobCores}:mem=${jobMemory} >> RUN/${scriptName}.sh
echo \#PBS -l walltime=${jobLength} >> RUN/${scriptName}.sh
echo \#PBS -o OUT/${scriptName}/ >> RUN/${scriptName}.sh
echo \#PBS -e ERR/${scriptName}/ >> RUN/${scriptName}.sh
echo \#PBS -J 1-${sampleNumber} >> RUN/${scriptName}.sh

echo module load ${moduleSTAR} >> RUN/${scriptName}.sh

echo sampleNames=\(${sampleNames[*]}\) >> RUN/${scriptName}.sh

## with default paramaters, STAR will map a read more than once, if there is ambiguity
## --outFilterMultimapNmax 10 is the default.  If more than 10 reads map, it'll be flaggfed as unmapped NH:0
## what happens here, is that we only report one of the multimap in the SAM file, so that the number of lines,
## is equal to the number of reads. --outSAMmultNMax 1 
## this is used in combination with --outMultimapperOrder Random which will pick a radnom top scoring read from 
## a set of reads with the same top score
## defaults output will not retained unmapped reads, pairs, this then screws up the bamtools statistic, therefor the following arguments is used: -outSAMunmapped Within KeepPairs

if [ ${fastqExtension} == "fq" ] || [ ${fastqExtension} == "fastq" ] ; then
	
	echo STAR --runThreadN 1 --runMode alignReads --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\.multimap\. --genomeDir ${myDir}/../resources/ --readFilesIn ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_1.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_2.${fastqExtension} >> RUN/${scriptName}.sh
#	echo STAR --outSAMunmapped Within KeepPairs --outMultimapperOrder Random --outSAMmultNmax 1 --runThreadN 1 --runMode alignReads --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\.onemap\. --genomeDir ${myDir}/../resources/ --readFilesIn ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_1.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_2.${fastqExtension} >> RUN/${scriptName}.sh
else

	## run star on trimmed data

	echo STAR --readFilesCommand zcat --runThreadN 1 --runMode alignReads --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\.multimap\. --genomeDir ${myDir}/../resources/ --readFilesIn ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_1.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_2.${fastqExtension} >> RUN/${scriptName}.sh
#	echo STAR --outSAMunmapped Within KeepPairs --outMultimapperOrder Random --outSAMmultNmax 1 --readFilesCommand zcat --runThreadN 1 --runMode alignReads --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${myDir}/../output/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}\.onemap\. --genomeDir ${myDir}/../resources/ --readFilesIn ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_1.${fastqExtension} ${myDir}/../tmp/\${sampleNames[\$PBS_ARRAY_INDEX - 1]}.trimmed\_2.${fastqExtension} >> RUN/${scriptName}.sh
	
fi

chmod u+x RUN/${scriptName}.sh 
qsub RUN/${scriptName}.sh


