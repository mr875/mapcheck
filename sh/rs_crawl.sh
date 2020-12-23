#!/bin/bash
# usage: rs_crawl.sh <rsid> 
baseurl='https://www.ncbi.nlm.nih.gov/snp/?term='
altbase='https://www.ncbi.nlm.nih.gov/snp/'
rsid=$1
url=$baseurl$rsid
out='out_rs_crawl'
wget "$url" -q -O $out 
grep "<title>Search: $rsid - NCBI</title>" $out > /dev/null
if [ $? -eq 0 ]; then
#    echo "$rsid not found"
    url2="$altbase$rsid"
    wget "$url2" -q -O $out 
    merges=$(sed -n '/<dt>Status<\/dt>/,/<\/dd>/p' $out | paste -s -d ' ' | sed 's/ \+/ /g' | sed 's/^.\+\(rs[0-9]\+\) \(was[^<]\+\)[^>]\+>\(rs[0-9]\+\).*$/\1 \3/')
    echo "merge: $merges"
    exit
fi
chrpos38=$(grep "</a></span></dd><dt>Chromosome: </dt><dd>" $out | sed 's/^.\+<dd>//' | uniq | paste -s -d ' ')
chrpos37=$(grep '(GRCh38)<br />' $out | sed 's/.\+\/>//' | uniq | paste -s -d ' ')
merges=$(egrep -o "rs[0-9]+ has merged into .+rs[0-9]+\">" out_rs_crawl | sed 's/^\(rs[0-9]\+\).\+\(rs[0-9]\+\).*$/\1 \2/' | paste -s -d ',')
echo -e "${rsid}\tb38=$chrpos38\tb37=$chrpos37\tMerges: $merges"
rm $out
