
# requires NCBI datasets CLI & python(v.3+) to be present in PATH

SPEC_LIST="bacteria_species.txt" # adjust or replace to include more species
OUTDIR="ref_genomes" # please adjust to a location fitting you


mkdir -p $OUTDIR/tmp
mkdir -p $OUTDIR/fnas

while read -r species; do
    species_noGap=$(echo $species | sed 's/ /_/')
    echo $species_noGap
    datasets download genome taxon "$species" --assembly-source 'RefSeq' --reference --filename $OUTDIR/tmp/ncbi_dataset.zip
    unzip $OUTDIR/tmp/ncbi_dataset.zip -d $OUTDIR/tmp
    mv $OUTDIR/tmp/ncbi_dataset/data/GCF*/*.fna $OUTDIR/fnas/
    rm -rf $OUTDIR/tmp/*
done < $SPEC_LIST

python createRefList.py --g $OUTDIR/fnas --o GARI_refList.tsv
