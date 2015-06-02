#!/bin/bash
# Converts pairwise IDR peak overlap output to narrowPeak

if [[ "$#" -lt 1 ]]
    then
    echo 'Converts pairwise IDR peak overlap output to narrowPeak' 1>&2
    echo "USAGE: $(basename $0) [idrOverlapFile] [oDir]" 1>&2
    echo '[idrOverlapFile]: overlap output file from pairwise IDR analysis' 1>&2
    echo '[oDir]: output directory' 1>&2
    exit 1
fi

# overlap file
ovFile=$1
if [[ ! -e ${ovFile} ]]
    then
    echo "ERROR:${ovFile} does not exist" 1>&2
    exit 1
fi

# Output directory
oDir=$(dirname ${ovFile})
[[ $# -gt 1 ]] && oDir=$2
if [[ ! -d ${oDir} ]]
    then
    mkdir ${oDir}
fi
oDir=$(echo ${oDir} | sed -r 's:/$::g')

# Create output file
oFile="${oDir}/$(basename ${ovFile} .gz).npk.gz"
if grep -q -E '\.gz$' ${ovFile}
    then
    zcat ${ovFile} | sed 1d | sed -r 's/"//g' | sort -k11g,11g | awk '{if ($3 <=$7) st=$3 ; else st=$7 ; if ($4 >= $8) sto=$4 ; else sto=$8 ; printf "%s\t%d\t%d\t%d\t%s\t.\t%s\t%f\t%f\n",$2,st,sto,NR,$5,$9,-log($10)/log(10),-log($11)/log(10)}' | gzip -c > ${oFile}
else
    sed 1d ${ovFile} | sed -r 's/"//g' | sort -k11g,11g | awk '{if ($3 <=$7) st=$3 ; else st=$7 ; if ($4 >= $8) sto=$4 ; else sto=$8 ; printf "%s\t%d\t%d\t%d\t%s\t.\t%s\t%f\t%f\n",$2,st,sto,NR,$5,$9,-log($10)/log(10),-log($11)/log(10)}' | gzip -c > ${oFile}
fi