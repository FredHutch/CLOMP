#!/bin/bash

#split -a 3 -d -l 100000 test.sam $basename

basename=$1

files=($basename*)
total=${#files[@]}
total=$((total-2))

echo "total"
echo $total

fixSplit() {
        inum=$(printf %03d $i)
        last_record=$(tail -1 $basename$inum | cut -f1)
        j=$((i + 1))
        j=$(printf %03d $j)
        echo "Last record in "$basename$inum" is "$last_record". Grepping in "$basename$j"."
        head -n 30000 $basename$j | grep $last_record >> $basename$inum
        echo "Removing "$last_record" from "$basename$j"."
        ex +:g/$last_record/d -cwq $basename$j
        echo "Done."
}

for i in $(seq 0 $total); do fixSplit & done

wait



# files=($basename*)
# total=${#files[@]}
# total=$((total-2))

# echo "total"
# echo $total

# fixSplit() {
#         inum=$(printf %03d $i)
#         echo $basename$inum
#         last_record=$(tail -1 $basename$inum | cut -f1)
#         echo $last_record
#         j=$((i + 1))
#         j=$(printf %03d $j)
#         echo "second file"
#         echo $basename$j
#         echo "startgrep"
#         head -n 30000 $basename$j | grep $last_record >> $basename$inum
#         echo "startsed"
#         ex +:g/$last_record/d -cwq $basename$j
#         echo "done"
# }

# for i in $(seq 0 $total); do fixSplit & done

# wait
