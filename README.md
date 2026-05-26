# ALS Drug Repositioning Pipeline

Takes a GWAS and returns a ranked drug table.
GWAS -> S-PrediXcan gene signature -> SignatureSearch (LINCS) -> approved-drug ranking.

Port of [`sarasaezALS/ALS-Drug-Repositioning`](https://github.com/sarasaezALS/ALS-Drug-Repositioning). Tested on OSC Ascend with the Van Rheenen 2021 ALS GWAS.

## Pipeline steps

```
config.yaml -> GWAS
   1. Harmonize       map to GTEx reference
   2. Format          add N, merge SE + frequency
   3. Impute          fill missing variants (220 jobs)
   4. Combine         merge imputed chunks
   5. S-PrediXcan     gene associations (signature tissue)
   6. Signature       spinal cord, FDR < 0.05, up/down
   7. SignatureSearch rank drugs by NCS
   8. Rank approved   filter to approved drugs, mean NCS
-> ranked approved-drug shortlist
```

Steps 1-5 run in the python env, 6-8 in the R env. This reproduces the drug-search trunk of the original analysis through the approved-drug ranking. S-MultiXcan (used only for a Miami plot in the original) is not included.

## Requirements

- A SLURM cluster
- conda
- ~80 GB disk for reference data and intermediates

## Get the code

Clone the repository and enter it. All commands below are run from the repository root.

```bash
git clone https://github.com/ggaitos/als_drug_repositioning_snakemake_pipeline.git
cd als_drug_repositioning_snakemake_pipeline
```

## Setup (once)

### Step 1 - Build the environments

```bash
conda create -p envs/snakemake -c conda-forge -c bioconda snakemake-minimal snakemake-executor-plugin-slurm -y
conda env create -p envs/imlabtools -f workflow/envs/imlabtools.yaml
conda env create -p envs/signaturesearch -f workflow/envs/signaturesearch.yaml
```

The R environment is a large Bioconductor solve; run it in tmux so a disconnect doesn't kill it.

### Step 2 - Get the reference data

From Zenodo record 3657902. If wget gives a 403, get the link from the API:

```bash
curl -s "https://zenodo.org/api/records/3657902" | python3 -c "import sys,json; d=json.load(sys.stdin); [print(f['links'].get('self')) for f in d['files']]"
cd resources && tar -xvf sample_data.tar && cd ..
```

This populates `resources/data/` with the coordinate map, 1000G panel, LD blocks, and MASHR models.

### Step 3 - Build the GENCODE table

```bash
mkdir -p resources/annotation && cd resources/annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
zcat gencode.v41.annotation.gtf.gz \
 | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' \
 | sed 's/gene_id "//; s/gene_type "//; s/gene_name "//; s/"//g' \
 | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' \
 | sed "1i\\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" \
 > gencode.v41_gene_annotation_table.txt
cd ../..
```

### Step 4 - Clone the upstream scripts

```bash
mkdir -p workflow/external && cd workflow/external
git clone https://github.com/hakyimlab/MetaXcan.git
git clone https://github.com/hakyimlab/summary-gwas-imputation.git
cd ../..
```

## Running your GWAS

### Step 5 - Add your GWAS and check it

Put your file in `resources/gwas/`, then check its header and build:

```bash
zcat resources/gwas/my_file.txt.gz | head -1
```

The pipeline expects hg19 with rsIDs. If your file is hg38, point `reference.coordinate_map` at `map_snp150_hg38.txt.gz`. If it has no rsID column, the quick harmonize step won't work as-is.

### Step 6 - Provide a drug list

The final step ranks only listed drugs. Put a plain text file at `resources/drugbank/approved_drugs.txt`, one drug name per line (lowercase). The ALS run used an approved-drug list of ~1000 names.

### Step 7 - Edit `config/config.yaml`

Left side is the role the pipeline needs; right side is your file's column name. Comments show the values used for the ALS run.

```yaml
gwas:
  name: MyGWAS                          # output goes to results/MyGWAS/
  file: resources/gwas/my_file.txt.gz
  columns:
    snp: <rsid column>                  # ALS: rsid
    non_effect_allele: <column>         # ALS: other_allele
    effect_allele: <column>             # ALS: effect_allele
    beta: <column>                      # ALS: beta
    pvalue: <column>                    # ALS: p_value
  sample_size: <total N>                # ALS: 138086
  n_cases: <number of cases>            # ALS: 27205
signature:
  tissue: <tissue>                      # ALS: Brain_Spinal_cord_cervical_c-1
drugbank:
  approved_list: resources/drugbank/approved_drugs.txt
```

### Step 8 - Run (inside tmux)

```bash
envs/snakemake/bin/snakemake -n                                  # dry run
envs/snakemake/bin/snakemake --workflow-profile profiles/slurm   # run
```

A new `name` writes to its own results folder, so other runs are untouched. Harmonize and format need ~32 GB; set your account and partition in `profiles/slurm/config.yaml`. Run inside tmux so a dropped connection doesn't kill the controller. If a run is interrupted, clear the leftover lock with `snakemake --unlock` before relaunching - it resumes from the last completed step.

## Output

In `results/<name>/`:

- `spredixcan/` gene associations for the signature tissue
- `signature/ALS.SpinalCord.Signature` the up/down signature
- `drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv` full ranked LINCS table
- `drugs/approved_mean_NCS.csv` approved drugs only, one row per drug, ranked by mean NCS (most negative = strongest reversal)

Check the signature has a known disease gene on top (C9orf72 for ALS):

```bash
column -t results/<name>/signature/ALS.SpinalCord.Signature | head
```

View the top approved candidates (already sorted):

```bash
head -20 results/<name>/drugs/approved_mean_NCS.csv
```

For the full LINCS table, parse as CSV rather than using `sort` (drug names contain commas):

```python
import csv
rows = [r for r in csv.DictReader(open("results/<name>/drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv")) if r["NCS"]]
rows.sort(key=lambda r: float(r["NCS"]))
for r in rows[:20]:
    print(f"{r['pert'][:30]:30} {r['cell']:6} {float(r['NCS']):8.3f}  {r['MOAss'][:40]}")
```

## Notes

- The approved-drug filter uses a plain name list, which substitutes for the original analysis's DrugBank join (the DrugBank tables aren't publicly redistributable). Matching is by exact lowercased name.
- Imputation parameters follow the MetaXcan CAD tutorial; the original ALS swarm file wasn't published, so they aren't confirmed against it.
- Python versions are pinned; R/Bioconductor is looser for portability, which can shift lower drug ranks but not top hits.
- Some output filenames keep `VR.ALS` from the original. A new GWAS still runs fine.
- Validated on Van Rheenen 2021: the spinal cord signature is led by C9orf72, and the approved-drug shortlist surfaces furosemide (validated experimentally in the original study) and edaravone (an approved ALS drug) among the top reversers.
