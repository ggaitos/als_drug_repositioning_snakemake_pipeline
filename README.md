Steps 1–5 run in the python env, 6–7 in the R env. This is the drug-search trunk of the original analysis; S-MultiXcan (Miami plot only) and the DrugBank approved-only filter are not included.

## Requirements

- A SLURM cluster
- conda
- ~80 GB disk for reference data and intermediates

## How to run it

### Step 1 — Build the environments (once)

```bash
conda create -p envs/snakemake -c conda-forge -c bioconda snakemake-minimal snakemake-executor-plugin-slurm -y
conda env create -p envs/imlabtools -f workflow/envs/imlabtools.yaml
conda env create -p envs/signaturesearch -f workflow/envs/signaturesearch.yaml
```

### Step 2 — Get the reference data (once)

From Zenodo record 3657902. If wget gives a 403, get the link from the API:

```bash
curl -s "https://zenodo.org/api/records/3657902" | python3 -c "import sys,json; d=json.load(sys.stdin); [print(f['links'].get('self')) for f in d['files']]"
cd resources && tar -xvf sample_data.tar
```

### Step 3 — Build the GENCODE table (once)

```bash
cd resources/annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation.gtf.gz
zcat gencode.v41.annotation.gtf.gz \
 | awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' \
 | sed 's/gene_id "//; s/gene_type "//; s/gene_name "//; s/"//g' \
 | awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' \
 | sed "1i\\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" \
 > gencode.v41_gene_annotation_table.txt
```

### Step 4 — Clone the upstream scripts (once)

```bash
cd workflow/external
git clone https://github.com/hakyimlab/MetaXcan.git
git clone https://github.com/hakyimlab/summary-gwas-imputation.git
```

### Step 5 — Add your GWAS and check it

Put your file in `resources/gwas/`, then check its header and build:

```bash
zcat resources/gwas/my_file.txt.gz | head -1
```

The pipeline expects hg19 with rsIDs. If your file is hg38, point `reference.coordinate_map` at `map_snp150_hg38.txt.gz`. If it has no rsID column, the quick harmonize step won't work as-is.

### Step 6 — Edit `config/config.yaml`

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
```

### Step 7 — Run (inside tmux)

```bash
envs/snakemake/bin/snakemake -n                                  # dry run
envs/snakemake/bin/snakemake --workflow-profile profiles/slurm   # run
```

A new `name` writes to its own results folder, so other runs are untouched. Harmonize and format need ~32 GB; set your account and partition in `profiles/slurm/config.yaml`.

## Output

In `results/<name>/`:

- `spredixcan/` per-tissue gene associations
- `signature/ALS.SpinalCord.Signature` the up/down signature
- `drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv` ranked drugs (most negative NCS = strongest reversal)

Check the signature has a known disease gene on top (C9orf72 for ALS):

```bash
column -t results/<name>/signature/ALS.SpinalCord.Signature | head
```

Sort the drug table as a CSV, not with `sort` (drug names contain commas):

```python
import csv
rows = [r for r in csv.DictReader(open("results/<name>/drugs/VR.ALS.lincsmethod.lincsDS.SpinalCordc1FDR.csv")) if r["NCS"]]
rows.sort(key=lambda r: float(r["NCS"]))
for r in rows[:20]:
    print(f"{r['pert'][:30]:30} {r['cell']:6} {float(r['NCS']):8.3f}  {r['MOAss'][:40]}")
```

## Notes

- Imputation parameters follow the MetaXcan CAD tutorial; the original ALS swarm file wasn't published, so they aren't confirmed against it.
- Python versions are pinned; R/Bioconductor is looser for portability, which can shift lower drug ranks but not top hits.
- Some output filenames keep `VR.ALS` from the original. A new GWAS still runs fine.
- Validated on Van Rheenen 2021: signature led by C9orf72, furosemide a top reverser.
