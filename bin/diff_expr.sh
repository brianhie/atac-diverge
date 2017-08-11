
cat $1 | \
    grep Yes | \
    grep ENSG | \
    cut -f $2 \
        > diff_expr.tmp

python bin/diff_expr.py \
       /godot/geuvadis/expression_analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz \
       conf/geuvadis_afr_eur.json \
       diff_expr.tmp | \
    sort -k4,4g -r | \
    head -n10
