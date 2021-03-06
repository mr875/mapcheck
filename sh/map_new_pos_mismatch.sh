#!/bin/bash
new_pos_mismatch=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
myname=${myname/map_/}
outfile="out_${myname%.*}_rs.txt"
while IFS=$'\t' read -r nothing id tabloc dblocs
do
    ori=${id#*orig }
    ori=${ori%)}
    id=${id%% *}
    if [[ $id =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $id)
        echo -e "$id\ttab=$tabloc\tdb=$dblocs\tbrorig=$ori\tlookup=$idlookup" | tee -a $outfile
    else
        echo -e "$id\ttab=$tabloc\tdb=$dblocs\tbrorig=$ori" | tee -a $outfile
    fi
done < $new_pos_mismatch
