#!/bin/bash
# $1: file containing arguments for batch-consistency-analysis.r
# File has 3 tab-delimited fields
# [peakFile1]\t[peakFile2]\t[pairOutFilePrefix]
# Rscript batch-consistency-analysis.r [peakfile1] [peakfile2] -1 [outfile.prefix] 0 F signal.value

if [[ "$#" -lt 1 ]]
then
	echo 'submit.idrpair.lsf.sh:'  1>&2
	echo 'Submits jobs to run pairwise IDR code'  1>&2
	echo 'USAGE:'  1>&2
	echo 'submit.idrpair.lsf.sh <idrPairArgFile>'  1>&2
	echo '<idrPairArgFile>: File has 3 tab-delimited fields' 1>&2
	echo '					[peakFile1]\t[peakFile2]\t[pairOutFilePrefix]' 1>&2
	exit 1
fi

rpath=`which Rscript`

JOBGROUPID="/idrPair${RANDOM}"

counter=1

while read inputline
do

  [[ $(( counter % 50 )) -eq 0 ]] && sleep 10s

  pf1="$(echo $inputline | awk '{print $1}')" # extract peak file 1
  pf2="$(echo $inputline | awk '{print $2}')" # extract peak file 2
  ofPrefix="$(echo $inputline | awk '{print $3}')" # extract outfile.prefix

  if [[ ! -e ${pf1} || ! -e ${pf2} || ! -d $( dirname ${ofPrefix} ) ]]
      then
      echo "Some file is not found for $( dirname ${ofPrefix} ): $(basename ${pf1}) $(basename ${pf2})"
      continue
  fi

  logfile="${ofPrefix}.log"
  randseed="${RANDOM}${RANDOM}"
  
  # If file exists then skip
  if [[ -e "${ofPrefix}-npeaks-aboveIDR.txt" ]]
      then
      continue
  fi

  # Create submit script
  submitScript="temp_${randseed}.sh"
  echo '#!/bin/bash' > ${submitScript}

  if echo ${pf1} | grep -q -E '\.gz$'
  then
    pf1gz=1
    pf1stub="$( basename ${pf1} | sed -r -e 's:\.gz$::g' )" # remove .gz and remove the directory name
    peakfile1="${TMP}/idr_${randseed}/${pf1stub}_${randseed}"
    echo "[[ ! -d ${TMP}/idr_${randseed} ]] && mkdir ${TMP}/idr_${randseed}" >> ${submitScript}
    echo "zcat ${pf1} > ${peakfile1}" >> ${submitScript}
  else
    pf1gz=0
    peakfile1="${pf1}"
  fi

  if echo ${pf2} | grep -q -E '\.gz$'
  then
    pf2gz=1
    pf2stub="$( basename ${pf2} | sed -r -e 's:\.gz$::g' )" # remove .gz and remove the directory name
    peakfile2="${TMP}/idr_${randseed}/${pf2stub}_${randseed}"
    echo "[[ ! -d ${TMP}/idr_${randseed} ]] && mkdir ${TMP}/idr_${randseed}" >> ${submitScript}
    echo "zcat ${pf2} > ${peakfile2}" >> ${submitScript}
  else
    pf2gz=0
    peakfile2="${pf2}"
  fi
  
  echo "${rpath} batch-consistency-analysis.r ${peakfile1} ${peakfile2} 500 ${ofPrefix} 0 F p.value" >> "${submitScript}"

  if [[ ${pf1gz} -eq 1 || ${pf2gz} -eq 1 ]]
      then
      echo "[[ -d ${TMP}/idr_${randseed} ]] && rm -rf ${TMP}/idr_${randseed}" >> ${submitScript}
  fi

  bsub -q research-rh6 -g ${JOBGROUPID} -W 48:00 -o ${logfile} -e ${logfile} < ${submitScript}
  (( counter = counter + 1 ))
  # bsub -g {JOBGROUPID} -W 48:00 -M 4096 -R "rusage[mem=4096]" -o ${logfile} -e ${logfile} < "${submitScript}"
  rm "${submitScript}"

done < $1

exit 0
