#!/bin/bash

paramFile="param.txt"
paramDir="/../resources/"

myDir=$(pwd)


# get parameter name from argument
paramName="$1"

if [ -z "${paramName}" ]; then
	exit
fi


# grep for parameter name in param file
paramLine=$(grep -P ${paramName}'\t' ${myDir}${paramDir}${paramFile})

if [ -n "${paramLine}" ]; then
	fieldValue=$(echo ${paramLine} | sed -e "s/${paramName}\s*//")
	echo $fieldValue
else
	exit
fi

