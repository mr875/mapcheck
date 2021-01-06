#!/bin/bash
new_alt_rs=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
myname=`basename "$0"`
outfile="out_${myname%.*}.txt"
while IFS=$'\t' read -r newalt line2 curr
do
    new=$(echo $newalt | egrep -o 'rs[0-9]+$')
    cur=$(echo $curr | egrep -o 'rs[0-9]+$')
    nlookup=$($crawl $new)
    clookup=$($crawl $cur)
    echo -e "New ($new):$nlookup\tCurrent ($cur):$clookup" | tee -a $outfile
done < $new_alt_rs
