#!/bin/bash
new_pos_mismatch=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
while IFS=$'\t' read -r nothing id tabloc dblocs
do
    id=${id%% *}
    dbloc=${dblocs%/*}
    echo $id $tabloc $dbloc
    if [[ $id =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $id)
        echo -e "$id\ttab=$tabloc\tdb=$dbloc\tlookup=$idlookup"
    else
        echo -e "$id\ttab=$tabloc\tdb=$dbloc"
    fi
#    cur=$(echo $curr | egrep -o 'rs[0-9]+$')
#    nlookup=$($crawl $new)
#    clookup=$($crawl $cur)
#    echo -e "New ($new):$nlookup\tCurrent ($cur):$clookup" 
done < $new_pos_mismatch
