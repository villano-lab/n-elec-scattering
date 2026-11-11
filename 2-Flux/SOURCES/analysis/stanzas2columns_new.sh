#! /bin/sh

file=$1

cat $file |grep -E 'Total Neutron Spectrum|Grand|Total S.F.|Multigroup' -A13 |awk 'BEGIN{a=1}{if($1+0==$1){for(i=1; i<NF; i++){printf("%1.4e  ",$i)}}else{printf("\n")}}' > test.txt

awk '{ for (i=1; i<=NF; i++) RtoC[i]= (i in RtoC?RtoC[i] OFS :"") $i; } END{ for (i=1; i<=NF; i++) print RtoC[i] }' test.txt
