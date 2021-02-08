#!/bin/bash
rs_byalt=$1
inpname=$(basename $rs_byalt)
inpname=${inpname%%_*}
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
myname=${myname/map_/}
outfile="$(date +%s)_${inpname}_out_${myname%.*}_rs.txt"
while IFS=$'\t' read -r newrs linkid mid
do
#    ori=${id#*orig };   ori=${ori%)}
    newrs=${newrs##* }
    linkid=${linkid##* }
    mid=${mid##* }
    idlookup=$($crawl $newrs)
    echo -e "map table rs= $newrs\tlink id= $linkid\tomics mid= $mid\tlookup= $idlookup" | tee -a $outfile
done < $rs_byalt
