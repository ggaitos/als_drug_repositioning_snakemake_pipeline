# ALS Drug Repositioning Pipeline

Snakemake pipeline that takes GWAS summary stats and returns a ranked drug table.
GWAS -> S-PrediXcan gene signature -> SignatureSearch (LINCS).

Port of `sarasaezALS/ALS-Drug-Repositioning`. Tested on OSC Ascend with the
Van Rheenen 2021 ALS GWAS.

## Steps

1. Harmonize GWAS to the GTEx reference
2. Format the harmonized file
3. Impute summary stats (220 jobs: 22 chr x 10 batches)
4. Combine imputed chunks
5. S-PrediXcan per tissue (whole blood + 13 brain)
6. Build signature from spinal cord (FDR < 0.05, up/down)
7. SignatureSearch vs LINCS, ranked by NCS

This covers the drug-search trunk of the original analysis. Two parts of the
notebook are outside this pipeline: S-MultiXcan (used only for a Miami plot) and
the final DrugBank approved-only filter (needs files that can't be redistributed).

## Requirements

- SLURM cluster (set for OSC Ascend, account PDE0075)
- conda
- ~80 GB disk

## Setup (once)

Build the environments on scratch:

    conda create -p envs/snakemake -c conda-forge -c bioconda snakemake-minimal snakemake-executor-plugin-slurm -y
    conda env create -p envs/imlabtools -f workflow/envs/imlabtools.yaml
    conda env create -p envs/signaturesearch -f workflow/envs/signaturesearch.yaml

Get the reference data (Zenodo record 3657902). If wget gives a 403, get the link from the API:

    curl -s "https://zenodo.org/api/records/3657902" | python3 -c "import sys,json; d=json.load(sys.stdin); [print(f['links'].get('self')) for f in d['files']]"
    cd resources && tar -xvf sample_data.tar

Build the GENCODE table:

    cd resources/annotation
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
    zcat gencode.v41.annotation.gtf.gz \
     | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' \
     | sed 's/gene_id "//; s/gene_type "//; s/gene_name "//; s/"//g' \
     | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' \
     | sed "1i\\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" \
     > gencode.v41_gene_annotation_table.txt

Clone the upstream scripts:

    cd workflow/external
    git clone https://github.com/hakyimlab/MetaXcan.git
    git clone https://github.com/hakyimlab/summary-gwas-imputation.git

## Configure

Edit `config/config.yaml` only. The rules don't change.

Check your GWAS header and build first:

    zcat resources/gwas/your_file.txt.gz | head -1

The pipeline assumes hg19 input with rsIDs. If your GWAS is hg38, point
`reference.coordinate_map` at `map_snp150_hg38.txt.gz` instead. If it has no rsID
column, the quick harmonize step won't work as-is.

Set in the config: `gwas.name`, `gwas.file`, the `gwas.columns` mapping,
`gwas.sample_size`, `gwas.n_cases`, and `signature.tissue` (spinal cord is
ALS-specific).

## Run

    envs/snakemake/bin/snakemake -n                                  # dry run
    envs/snakemake/bin/snakemake --workflow-profile profiles/slurm   # run

Run inside tmux. A new `gwas.name` writes to its own results folder. Harmonize and
format need ~32 GB; set account and partition in `profiles/slurm/config.yaml`.

## Output

In `results/<gwas.name>/`:

- `spredixcan/` per-tissue gene associations
- `signature/ALS.SpinalCord.Signature` the up/down signature
- `drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv` ranked drugs

Most negative NCS means strongest predicted reversal of the signature.

Check the signature has a known disease gene on top (C9orf72 for ALS):

    column -t results/<gwas.name>/signature/ALS.SpinalCord.Signature | head

Sort the drug table as a CSV, not with `sort` (drug names contain commas):

    envs/imlabtools/bin/python - <<'PY'
    import csv
    rows = [r for r in csv.DictReader(open("results/<gwas.name>/drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv")) if r["NCS"]]
    rows.sort(key=lambda r: float(r["NCS"]))
    for r in rows[:20]:
        print(f"{r['pert'][:30]:30} {r['cell']:6} {float(r['NCS']):8.3f}  {r['MOAss'][:40]}")
    PY

## Notes

- Imputation parameters follow the MetaXcan CAD tutorial; the original ALS swarm
  file wasn't published, so they aren't confirmed against it.
- Python versions are pinned (pyarrow and pandas matter); R/Bioconductor is looser
  for portability, which can shift lower drug ranks but not the top hits.
- Some output filenames keep `VR.ALS` from the original. A new GWAS still runs fine.
- Validated on Van Rheenen 2021: signature led by C9orf72, furosemide a top reverser.
