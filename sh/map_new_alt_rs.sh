#!/bin/bash
new_alt_rs=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
#awk -v FS='\t' '{print $1,$3}' $new_alt_rs | sed 's/^.\+(uid): \(rs[0-9]\+\) .\+id: \(rs[0-9]\+\)/\1 \2/'
echo $crawl
while IFS=$'\t' read -r newalt line2 curr
do
    new=$(echo $newalt | egrep -o 'rs[0-9]+$')
    cur=$(echo $curr | egrep -o 'rs[0-9]+$')
    nlookup=$($crawl $new)
    clookup=$($crawl $cur)
    echo -e "New ($new):$nlookup\tCurrent ($cur):$clookup" 
done < $new_alt_rs
