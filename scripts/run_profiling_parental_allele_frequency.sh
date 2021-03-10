batch="Batch_4viable_top20_6";
mode="raw"

perl batch_parental_allele_frequency_in_tetrads_profiling.pl -s Master_Sample_Table.$batch.txt -genotype_dir $batch -m $mode -o $batch.parental_allele_frequency.$mode.txt

Rscript  --vanilla --slave plot_parental_allele_frequency_in_tetrads.R \
    --input $batch.parental_allele_frequency.$mode.txt \
    --output $batch.parental_allele_frequency.$mode.plot.pdf \
    --color_scheme ./../../data/color_scheme.4way.txt


