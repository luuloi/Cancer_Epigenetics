# Dataset Description — Osmotic Stress ATAC-seq Time Course

**Author:** Quynh Nhu Nguyen

---

## Sample Information

| Field | Details |
|---|---|
| **Organism** | *Saccharomyces cerevisiae* (baker's yeast) |
| **Strain** | BY4741 |
| **Growth stage** | Mid-log phase |
| **Growth medium** | YPD |
| **Published** | April 22, 2015 |
| **Source** | NCBI SRA |
| **Instrument** | Illumina NextSeq 500 |
| **Genome assembly** | sacCer3 (*S. cerevisiae*) |

---

## Data

**Raw data storage:** [Google Drive — Raw FASTQ files](https://drive.google.com/drive/u/0/folders/1ODFlrXsItKVkucM-Pk-xItjGYkzMkdDL)

| Sample | Replicate | R1 | R2 |
|---|---|---|---|
| OSMOTIC_STRESS_T0_PE | rep1 | SRR1822153_1.fastq.gz | SRR1822153_2.fastq.gz |
| OSMOTIC_STRESS_T0_PE | rep2 | SRR1822154_1.fastq.gz | SRR1822154_2.fastq.gz |
| OSMOTIC_STRESS_T15_PE | rep1 | SRR1822157_1.fastq.gz | SRR1822157_2.fastq.gz |
| OSMOTIC_STRESS_T15_PE | rep2 | SRR1822158_1.fastq.gz | SRR1822158_2.fastq.gz |

> **Note:** Files in `raw/` are downsampled to 100,000 reads per file for pipeline testing only. Full raw data is available at the Google Drive link above.

---

## Experimental Design

**Experiment:** ATAC-seq measuring changes in chromatin accessibility in response to **osmotic stress** — a condition where high external salt concentration causes water to leave the cell, leading to cell shrinkage and triggering a stress response.

| Condition | Time point | Treatment | Meaning |
|---|---|---|---|
| `OSMOTIC_STRESS_T0` | T = 0 min | No NaCl | Normal chromatin — baseline |
| `OSMOTIC_STRESS_T15` | T = 15 min | 0.6 M NaCl added | Chromatin after 15 min of osmotic stress |

Each condition has **2 biological replicates**.

---

## What is Osmotic Stress?

Osmotic stress occurs when the external environment has a higher salt concentration than the inside of the cell. Water flows out of the cell down the osmotic gradient, causing the cell to shrink.

```
Normal:
  [Cell]  ←→  [YPD medium]
  balanced pressure → healthy, turgid cell

After adding 0.6 M NaCl:
  [Cell]  →→→  [Salty medium]
  water exits the cell along the osmotic gradient
  → cell shrinks (plasmolysis)
```

### Cell response during the first 15 minutes

| Time | Response |
|---|---|
| 0–2 min | Water exits, cell shrinks, chromatin compacts |
| 2–10 min | Cell activates the **HOG pathway** (High Osmolarity Glycerol) — an emergency signaling cascade |
| 10–15 min | Stress response genes begin to be transcribed; nucleosomes are evicted from promoters to allow RNA polymerase access |

### Why T = 15 minutes?

This is the moment of **peak active response** — the cell has not yet fully adapted but has already begun opening chromatin at genes needed for survival. At T = 60 minutes the cell has fully adapted and chromatin returns nearly to baseline, making changes harder to detect.

---

## Biological Question

> When *S. cerevisiae* is subjected to osmotic stress (0.6 M NaCl), how do nucleosomes reposition? Which genomic regions become more accessible, and which close down? Which genes are activated first?

Comparing **T0 vs T15** reveals precisely which regions of the genome open up during the first 15 minutes of stress, identifying the earliest transcriptional response to high-salt conditions.

---

## Extraction Protocol

1. Harvest **5 million cells** at mid-log phase
2. Wash once in **Sorbitol buffer** (1.4 M Sorbitol, 40 mM HEPES-KOH pH 7.5, 0.5 mM MgCl2) + 10 mM DTT + 0.6 M NaCl *(T15 only)*
3. Incubate 5 min at 30°C, 300 rpm with **0.5 mg/mL zymolyase** to digest the cell wall
4. Wash once more in Sorbitol buffer (+ 0.6 M NaCl for T15)
5. Add **2.5 µL Nextera Transposase** in 47.5 µL 1× TD buffer; incubate 37°C, 300 rpm, 15 min
6. Construct sequencing libraries following **Buenrostro et al. (2013)**

---

## Data Processing (as described by original authors)

- Adapter trimming before alignment
- Alignment with **Bowtie2**, flag `-X 2000` (maximum insert size 2000 bp)
- Discard: duplicate fragments, reads mapping to chrM, MAPQ < 30, improperly paired reads
- Open chromatin peak calling with **MACS2**
- Nucleosome positioning within open chromatin regions with **NucleoATAC**

---

## References

### Dataset

Schep AN, Buenrostro JD, Cagney G, Garber M, Haber JE, Bhatt DL, **Greenleaf WJ**.
*Structured nucleosome fingerprints enable high-resolution mapping of chromatin architecture within regulatory regions.*
**Genome Research**, 2015.

| Field | Value |
|---|---|
| GEO Series | [GSE66386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66386) |
| GEO Title | High-resolution nucleosome positioning from ATAC-seq chromatin accessibility data |
| Accession | PRJNA276699; GEO: GSE66386 |
| Submission date | Feb 27, 2015 |
| Contact | William J. Greenleaf — wjg@stanford.edu |
| Institution | Stanford University, Department of Genetics |
| SRA accessions | SRR1822153, SRR1822154, SRR1822157, SRR1822158 |
| Genome build | sacCer3 |

