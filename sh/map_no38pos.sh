#!/bin/bash
no38pos=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
crawl="${DIR}/rs_crawl.sh"
while IFS=$'\t' read -r nothing id pos
do
    id=${id%% *}
    if [[ $id =~ ^rs[0-9] ]]; then
        idlookup=$($crawl $id)
        echo -e "$id\t$pos\tlookup:$idlookup"
    else
        echo -e "$id\t$pos"
    fi
done < $no38pos
