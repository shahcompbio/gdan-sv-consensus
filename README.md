# SV consensus pipeline
Created to call consensus SVs from the SV calls of GDAN DLBCL groups.

# Usage
1. Fix `config.yaml`

- `samples_file`, `metadata` has to be generated
- `ref` should be relinked to your local references

2. Install prerequisites

Haven't set up a `setup.py` yet, you could just type in
```bash
pip install tabix numpy pandas pybedtools wgs_analysis
```

3. Fix global variables in `Snakefile`

You may be interested in fixing `CHROMS`, `SOURCES`, `SAMPLES` to your liking.

4. Run

Now the fun part! A snakemake dry-run is highly recommended.
```
bash run_snakemake.sh
```

# BEDPE output explanation
The output BEDPE files created by this pipeline has the following columns:
- `#chrom1`: chromosome of breakpoints 1
- `start1`, `end1`: position of breakpoint 1, in 0-based semi-noninclusive range (i.e. if a 1-based UCSC breakpoint is `chr1:1`, the breakpoint coordinate in the BEDPE is `chr1 0 1`
- `chrom2`, `start2`, `end2`: coordinates of breakpoint 2
- `name`: sources of the SV separated by `__`
- `score`: number of the sources of the SV
- `strand1`, `strand2`: orientation of breakpoint 1 and 2, respectively
- `type`: type of the SV, in {'DEL', 'DUP', 'INV', 'TRA'}
- `gene1`, `gene2`: gene annotation of breakpoint 1 and 2

## Gene annotation of a breakpoint
1. For `DEL` and `DUP`
- For breakpoint 1, closest downstream gene that is encompassed by the deletion/duplication (i.e. has overlapping coordinates with the breakpoint pair interval) is annotated
- For breakpoint 2, closest upstream gene that is encompassed by the deletion/duplication is annotated
2. For `INV` and `TRA`
- If a breakpoint orientation is upstream (i.e. '+'), the closest upstream gene to the breakpoint is annotated
- If a breakpoint orientation is downstream (i.e. '-'), the closest downstream gene to the breakpoint is annotated
