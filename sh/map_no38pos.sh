#!/bin/bash
no38pos=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
myname=${myname/map/out}
outfile="${myname%.*}_rs.txt"
while IFS=$'\t' read -r nothing id pos
do
    id=${id%% *}
    if [[ $id =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $id)
        echo -e "$id\t$pos\tlookup:$idlookup" | tee -a $outfile
    else
        echo -e "$id\t$pos" | tee -a $outfile
    fi
done < $no38pos
