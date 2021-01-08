#!/bin/bash
mismatflag=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
outfile="out_${myname%.*}_rs.txt"
while IFS=$'\t' read -r noth id pos
do
    pid=${id%% *}
    oid=${id##* }
    oid=${oid::-1}
    if [[ $pid =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $pid)
        echo -e "$pid\torig:$oid\tlookup:$idlookup" | tee -a $outfile
    else
        echo -e "$pid\torig:$oid" | tee -a $outfile
    fi

done < $mismatflag
