# ALS Drug Repositioning Pipeline

Snakemake pipeline: GWAS summary stats -> S-PrediXcan gene signature -> SignatureSearch (LINCS) drug ranking.

Port of `sarasaezALS/ALS-Drug-Repositioning`. Built and tested on OSC Ascend with the Van Rheenen 2021 ALS GWAS.

## Steps

1. Harmonize GWAS to GTEx reference (MetaXcan M03_betas.py)
2. Format harmonized file (add N, merge back SE + frequency)
3. Impute summary stats (220 jobs: 22 chr x 10 batches)
4. Combine imputed chunks
5. S-PrediXcan per tissue (14: whole blood + 13 brain)
6. Build signature from spinal cord (FDR < 0.05, up/down split)
7. SignatureSearch vs LINCS -> ranked drugs

## Requirements

- SLURM cluster (set up for OSC Ascend, account PDE0075)
- conda
- ~80 GB disk for reference data + intermediates

## Setup (one time)

**Environments** (build on scratch with `-p` due to home quota):

    conda create -p envs/snakemake -c conda-forge -c bioconda snakemake-minimal snakemake-executor-plugin-slurm -y
    conda env create -p envs/imlabtools -f workflow/envs/imlabtools.yaml
    conda env create -p envs/signaturesearch -f workflow/envs/signaturesearch.yaml

The R env is a big Bioconductor solve; run it in tmux. (YAMLs are derived from a working env, not rebuilt from scratch elsewhere — a fresh solve may need minor adjustment.)

**Reference data** (not auto-downloaded). Get sample_data.tar from Zenodo record 3657902. If wget 403s, get the real link from the API:

    curl -s "https://zenodo.org/api/records/3657902" | python3 -c "import sys,json; d=json.load(sys.stdin); [print(f['links'].get('self')) for f in d['files']]"

Then:

    cd resources && tar -xvf sample_data.tar

Gives `resources/data/` with the coordinate map, 1000G panel (per-chr parquet), LD blocks, and MASHR models for 49 tissues.

**GENCODE table** for the signature step:

    cd resources/annotation
    wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
    zcat gencode.v41.annotation.gtf.gz \
     | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' \
     | sed 's/gene_id "//; s/gene_type "//; s/gene_name "//; s/"//g' \
     | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' \
     | sed "1i\\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" \
     > gencode.v41_gene_annotation_table.txt

**Upstream scripts:**

    cd workflow/external
    git clone https://github.com/hakyimlab/MetaXcan.git
    git clone https://github.com/hakyimlab/summary-gwas-imputation.git

## Configure for your GWAS

Edit only `config/config.yaml` — the rules never change.

**First, check your GWAS file.** Look at the header and confirm the genome build:

    zcat resources/gwas/your_file.txt.gz | head -1

This pipeline assumes **hg19** input with **rsIDs** (harmonize maps rsIDs via map_snp150_hg19). Two checks:
- If your GWAS is hg38, point `reference.coordinate_map` at `map_snp150_hg38.txt.gz` (included in the sample data) instead.
- If your file has no rsID column (only chr:pos), the quick harmonize step won't work as-is and needs the fuller gwas_parsing.py approach — not covered here.

**Then edit config/config.yaml:**

    gwas:
      name: MyGWAS                       # sets output folder: results/MyGWAS/
      file: resources/gwas/my_file.txt.gz
      columns:                           # map YOUR headers to these roles
        snp: <rsid column>
        non_effect_allele: <column>
        effect_allele: <column>
        beta: <column>
        pvalue: <column>
      sample_size: <N>
      n_cases: <N>

    signature:
      tissue: <relevant tissue>          # spinal cord is ALS-specific; change per disease

Imputation params (window, frequency_filter, regularization, sub_batches) are MetaXcan CAD-tutorial defaults; the original ALS swarm file wasn't published.

## Run

    envs/snakemake/bin/snakemake -n                                  # dry run
    envs/snakemake/bin/snakemake --workflow-profile profiles/slurm   # run on SLURM

Run inside tmux. Changing `gwas.name` writes to a new `results/` folder, so previous runs are untouched and a new GWAS runs from scratch. SLURM account/partition/memory live in `profiles/slurm/config.yaml`; harmonize and format need ~32 GB.

## Output

In `results/<gwas.name>/`:
- `spredixcan/` — per-tissue gene associations
- `signature/ALS.SpinalCord.Signature` — up/down gene signature
- `drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv` — ranked drugs (most negative NCS = strongest reversal)

Some output filenames contain `VR.ALS` literally (kept from the original workflow). A new GWAS still runs correctly; the files just keep those names.

## Check the output

Signature — expect a known disease gene near the top (C9orf72 for ALS):

    column -t results/<gwas.name>/signature/ALS.SpinalCord.Signature | head

Drug table — parse as CSV, don't `sort` (drug names have commas):

    envs/imlabtools/bin/python - <<'PY'
    import csv
    rows = [r for r in csv.DictReader(open("results/<gwas.name>/drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv")) if r["NCS"]]
    rows.sort(key=lambda r: float(r["NCS"]))
    for r in rows[:20]:
        print(f"{r['pert'][:30]:30} {r['cell']:6} {float(r['NCS']):8.3f}  {r['MOAss'][:40]}")
    PY

## Known limitations

- Final DrugBank approved-only ranking not included (needs non-redistributable DB.drug_groups.tab, information.db.rda). Pipeline outputs the full LINCS table.
- S-MultiXcan not run (fed a Miami plot in the original, not the drug list).
- Python pins are exact (pyarrow/pandas matter); R/Bioconductor left looser for portability — a fresh build may shift lower drug ranks, not the top biology.

## Note

Validated on Van Rheenen 2021. Signature led by C9orf72; furosemide appears with strongly negative NCS, matching the original study.
