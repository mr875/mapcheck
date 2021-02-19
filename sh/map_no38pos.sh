#!/bin/bash
no38pos=$1
inpname=$(basename $no38pos)
inpname=${inpname%%_*}
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
myname=${myname/map_/}
#outfile="${myname%.*}_rs.txt"
outfile="$(date +%s)_${inpname}_out_${myname%.*}_rs.txt"
while IFS=$'\t' read -r nothing id pos
do
    ori=${id##* }
    ori=${ori%)}
    id=${id%% *}
    if [[ $id =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $id)
        echo -e "$id\t$pos\tbrorig=$ori\tlookup:$idlookup" | tee -a $outfile
    else
        echo -e "$id\t$pos\tbrorig=$ori" | tee -a $outfile
    fi
done < $no38pos
