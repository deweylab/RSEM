#!/bin/bash
# $1: file containing arguments for batch-consistency-plot-merged3.r
# File has 4 tab-delimited fields
# [nPairs]\t[combOutFilePrefix]\t[pairOutFilePrefix,...]\t[combPeakFile]
# $2: idrThreshold (OPTIONAL: default is 0.1)
# $3: fdrThreshold (OPTIONAL: default is 0.7)
# Rscript batch-consistency-plot-merged2.r [npairs] [output.prefix] [input.file.prefix 1, 2, 3 ...] -1 [idr.level] signal.value [pooled.filename] T 0 F

if [[ "$#" -lt 1 ]]
then
	echo 'submit.idrmerge.lsf.sh:'  1>&2
	echo 'Submits jobs to run IDR code on pooled peak calls'  1>&2
	echo 'USAGE:'  1>&2
	echo 'submit.idrmerge.lsf.sh <idrMergeArgFile> <OPTIONAL:idrThresh> <OPTIONAL:fdrThresh>'  1>&2
	echo '<idrMergeArgFile>: File has 4 tab-delimited fields' 1>&2
	echo '					[nPairs]\t[combOutFilePrefix]\t[pairOutFilePrefix,...]\t[combPeakFile]' 1>&2
	echo '<idrThresh>: OPTIONAL: Default of 0.1' 1>&2
	echo '<fdrThresh>: OPTIONAL: Default of 0.7' 1>&2
	exit 1
fi

idrlevel='0.1'
if [[ "$#" -ge 2 ]]; then idrlevel=$2 ; fi     
fdrthresh='0.7'
if [[ "$#" -ge 3 ]]; then fdrthresh=$3 ; fi     
fdrthresh=$(echo ${fdrthresh} | awk '{print -log($1)/log(10)}') # convert fdrthreshold or logscale

rpath=`which Rscript`
# TEMPORARY DIRECTORY
#TMP='/scratch/users/akundaje/temp'

while read inputline
do

  npairs="$(echo $inputline | awk '{print $1}')" # extract npairs
  ofPrefix="$(echo $inputline | awk '{print $2}')" # extract output.prefix
  ifPrefix="$(echo $inputline | awk '{print $3}' | sed -r 's/,/ /g')" # extract input.prefixes
  combFname="$(echo $inputline | awk '{print $4}')" # extract merged peak file name
  echo "${ofPrefix}"
  logfile="${ofPrefix}.log"
  randseed="$RANDOM"

  # Create submitScript
  submitScript="tempMerge_${randseed}.sh"
  echo '#!/bin/bash' > "${submitScript}"  

  if [[ `echo ${combFname} | grep -E '\.gz$'` ]]
  then
    isgz='1'
    combStub="$(echo ${combFname} | sed -r -e 's:\.gz$::g' -e 's:^.*/::g')" # remove .gz and remove the directory name
    combPeakFile="${TMP}/${combStub}_${randseed}"
    # echo "gunzip -c ${combFname} > ${combPeakFile}" >> "${submitScript}"
  else
    isgz='0'
    combPeakFile="${combFname}"
  fi
 
  combPeakFile='random.txt'
  echo "${rpath} batch-consistency-plot-merged2.r ${npairs} ${ofPrefix} ${ifPrefix} -1 ${idrlevel} signal.value ${combPeakFile} T 0 F" >> "${submitScript}"

  if [[ "${isgz}" == '1' ]]
  then
    echo "rm -f ${combPeakFile}" >> "${submitScript}"
  fi
  
  chmod 755 "${submitScript}"
  bsub -W 24:00 -M 4096 -R "rusage[mem=4096]" -o ${logfile} -e ${logfile} < "${submitScript}"    
  rm "${submitScript}"

done < $1

exit 0
