#!/bin/bash
new_rs=$1
inpname=$(basename $new_rs)
inpname=${inpname%%_*}
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
myname=${myname/map_/}
outfile="$(date +%s)_${inpname}_out_${myname%.*}_rs.txt"
while IFS=$'\t' read -r newrs secid
do
#    ori=${id#*orig };   ori=${ori%)}
    newrs=${newrs##* }
    secid=${secid##* }
    idlookup=$($crawl $newrs)
    echo -e "map table rs= $newrs\tlink id= $secid\tlookup= $idlookup" | tee -a $outfile
done < $new_rs
