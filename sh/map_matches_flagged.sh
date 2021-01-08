#!/bin/bash
matflag=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
if [[ $matflag == *"alt_avail"*  ]]; then
    outfile="out_${myname%.*}_alt_avail_rs.txt"
else
    outfile="out_${myname%.*}_rs.txt"
fi    

while IFS=$'\t' read -r id pos comp
do
    pid=${id%% *}
    oid=${id%%)*}
    oid=${oid##* }
    if [[ $pid =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $pid)
        echo -e "$pid\torig:$oid\t$comp\tlookup:$idlookup" | tee -a $outfile
    else
        echo -e "$pid\torig:$oid\t$comp" | tee -a $outfile
    fi
done < $matflag
