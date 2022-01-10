#!/bin/bash

catch=''

pathIn="/data/stohyd/mHM_project/germany/basin/sub_${catch}/output/params/"
parameterFile="/data/stohyd/mHM_project/germany/basin/sub_${catch}/output/${catch}_parameter.sets"

cd ${pathIn}

# grep all files with optimized parameters
# print every second line
# pipe it into file 'tmp'
cat $(find ${pathIn} -name 'FinalParam.out') | sed -n '2~2p' > tmp
# number of parameters which are already in the parameterFile
let NoParOld=$(wc -l ${parameterFile} | cut -d ' ' -f 1)-1
let NoParNew=$(wc -l tmp | cut -d ' ' -f 1)

echo 'Starting loop'
# check if tmp2 file exists
if [ -f tmp2 ] ; then rm tmp2 ; fi
#for (( i=1 ; i < 4 + 1 ; i++ )) ; do
for (( i=1 ; i < ${NoParNew} + 1 ; i++ )) ; do
    let no=${NoParOld}+${i}
    echo 'Processing parameter ...' $no
    # print in the same format like parameterFile
    printf "%7s%13s%13s%15s%15s"\
            ${no} $(sed -n "${i}p" tmp | tr -s ' ' ' '| cut -d ' ' -f 3-6 ) >> tmp2
    printf "%12s"\
    $(sed -n "${i}p" tmp | tr -s ' ' ' '| cut -d ' ' -f 7- ) >> tmp2
    printf "\n" >> tmp2
done

cat ${parameterFile} tmp2 > tmp

mv ${parameterFile} ${parameterFile}'.old'

mv tmp ${parameterFile}

rm tmp2

exit 0
