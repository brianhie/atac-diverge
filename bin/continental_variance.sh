#!/usr/bin/bash

TYPE=reads

if [ ! -e target/pop_peak_"$TYPE".txt ]
then
    bin/main.sh
fi

echo `date`" | Finding peaks with continental divergence..."
python bin/continental_variance.py target/pop_peak_"$TYPE".txt \
       > target/continental_variance/"$TYPE".txt

echo `date`" | Separating divergence into African and European biased..."
awk '$5 > $6' target/continental_variance/"$TYPE".txt \
    > target/continental_variance/"$TYPE"_biased_afr.txt
awk '$5 < $6' target/continental_variance/"$TYPE".txt \
    > target/continental_variance/"$TYPE"_biased_eur.txt

echo `date`" | Mapping peaks to rsIDs..."
CONTINENTS=(afr eur)
for CONTINENT in ${CONTINENTS[@]}
do
    python bin/peak_to_rsid.py \
           data/snps_in_peaks.bed \
           target/continental_variance/"$TYPE"_biased_"$CONTINENT".txt \
           > target/continental_variance/"$TYPE"_biased_"$CONTINENT"_rsids.txt &
done
wait

echo `date`" | Done!"
