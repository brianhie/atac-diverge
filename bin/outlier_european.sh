#!/usr/bin/bash

OUTLIER_POP=yri
TYPE=reads

if [ ! -e target/pop_peak_"$TYPE".txt ]
then
    bin/main.sh
fi

#echo `date`" | Finding peaks with outlier population..."
#python bin/outlier_european.py target/pop_peak_"$TYPE".txt $OUTLIER_POP \
#       > target/outlier_european/"$TYPE".txt
#
#echo `date`" | Separating based on bias..."
#awk '$5 > $6' target/outlier_european/"$TYPE".txt \
#    > target/outlier_european/"$TYPE"_biased_"$OUTLIER_POP".txt
#awk '$5 < $6' target/outlier_european/"$TYPE".txt \
#    > target/outlier_european/"$TYPE"_biased_eur.txt

echo `date`" | Mapping peaks to TSSs..."
POPS=($OUTLIER_POP eur)
for POP in ${POPS[@]}
do
    (
        cat target/outlier_european/"$TYPE"_biased_"$POP".txt | \
            sed 's/\t/\./' | sed 's/\t/\./' | \
            sort -k1,1 \
                 > temp_"$TYPE"_"$POP"-outlier_european.txt.1

        cat target/pop_peak_"$TYPE".txt | \
            sed 's/\t/\./' | sed 's/\t/\./' | \
            sort -k1,1 | \
            join - temp_"$TYPE"_"$POP"-outlier_european.txt.1 | \
            sed 's/\./\t/' | sed 's/\./\t/' | sed 's/ /\t/g' \
                > temp_"$TYPE"_"$POP"-outlier_european.txt.2
        
        python bin/correlate_peak_expr.py \
               data/genes_unique.txt \
               temp_"$TYPE"_"$POP"-outlier_european.txt.2 \
               /godot/geuvadis/expression_analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz \
               conf/geuvadis_pops.json \
               > target/outlier_european/"$TYPE"_biased_"$POP"_genes.txt

        rm -rf temp_"$TYPE"_"$POP"-*
    ) &
done
wait

echo `date`" | Done!"
