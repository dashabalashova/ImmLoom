# ImmLoom

Pipeline for IG locus analysis using patchworkplot + graphs + MSA.

## Run everything

```bash
./runs/run_dataset-02-IGH.sh
```

## What it does

Pipeline steps:

1. Setup env + install deps  
2. Download + prepare dataset  
3. Run patchworkplot  
4. Fix sequence orientation  
5. Build and connect blocks  
6. Build graph components
7. Export FASTA  
8. Run Clustal Omega  
9. Build MSA figures  

## Outputs

Main results are in:

```text
results/dataset-02/IGH/
```
