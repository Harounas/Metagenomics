while IFS= read -r line
do
    # Split the line into array to get paired-end files
    files=($line)
    #s0="${files%_*}"
 #dir_path=$(dirname "$files") 
 s1=$(basename "$files")
s2=${s1%.*}
#echo ${dir_path}
echo $s2

#awk -v OFS=, '{print  $5}' ${s2}.tsv|sort|uniq -c|awk 'BEGIN {OFS="\t"} {print $2, $1}'>${s2}.1tsv 
cut -f4 OSL012_S34_kaiju.names.out | awk -v OFS="\t" '{print $0}'|sort|uniq -c|grep -v -i "NA;">${s2}.otsv
awk -v OFS="\t" '{ printf "%s %s %s %s %s %s\t%s\n",$2,$3,$4,$5,$6,$7,$1; }' ${s2}.otsv>${s2}.tsv
done < "summaryfiles.txt"


