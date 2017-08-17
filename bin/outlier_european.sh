#!/usr/bin/bash

OUTLIER_POP=yri
TYPE=reads

if [ ! -e target/pop_peak_"$TYPE".txt ]
then
    bin/main.sh
fi

echo `date`" | Finding peaks with outlier population..."
python bin/outlier_european.py target/pop_peak_"$TYPE".txt $OUTLIER_POP \
       > target/outlier_european/"$TYPE".txt

echo `date`" | Separating based on bias..."
awk '$5 > $6' target/outlier_european/"$TYPE".txt \
    > target/outlier_european/"$TYPE"_biased_"$OUTLIER_POP".txt
awk '$5 < $6' target/outlier_european/"$TYPE".txt \
    > target/outlier_european/"$TYPE"_biased_eur.txt

echo `date`" | Mapping peaks to rsIDs..."
POPS=($OUTLIER_POP eur)
for POP in ${POPS[@]}
do
    python bin/peak_to_rsid.py \
           depict/data/trityper_CEU_hg19/SNPMappings.txt \
           target/outlier_european/"$TYPE"_biased_"$POP".txt \
           > target/outlier_european/"$TYPE"_biased_"$POP"_rsids.txt &
done
wait

echo `date`" | Running DEPICT..."
for POP in ${POPS[@]}
do
    cut -f1 target/outlier_european/"$TYPE"_biased_"$POP"_rsids.txt \
        > depict/testfiles/"$TYPE"_biased_"$POP"_rsids.txt
    (
        cd depict;
        ./depict.py "$TYPE"_biased_"$POP" outlier_european
    ) &
done
wait

echo `date`" | Done!"
