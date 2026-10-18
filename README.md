# SOG1 in banana: phase 2 of my *Ma-SOG1* project

**Live report:** https://bymankirat.github.io/SOG1-Banana-Comparative-Genomics/

SOG1 (SUPPRESSOR OF GAMMA RESPONSE 1) is the plant counterpart of p53. After DNA damage the ATM and ATR kinases phosphorylate
SQ/TQ sites in its tail, and SOG1 then switches on DNA repair, stops the cell cycle or triggers cell death.

In 2023, during research training at the National Agri-Food Biotechnology Institute (NABI), Mohali, I analysed the banana gene
*Ma-SOG1* in silico and cloned it into a plant expression vector for overexpression in *Arabidopsis*. That lab work belongs
to the NABI group and is not included here. This repository is the computational next step, built only on public genomes
(17 plant genomes from Ensembl Plants release 63, from green algae to banana).

## What I found

| Question | Result |
|---|---|
| Which banana gene was cloned in phase 1? | **Ma03_g07280**, the one-to-one ortholog of Arabidopsis SOG1 (also the Ensembl Compara call) |
| How many SOG1-clade genes does banana have? | **4** of its 181 NAC genes: Ma03_g07280 (SOG1), Ma00_g02930 and Ma09_g01850 (NAC044-like, a paralog pair, Ks 0.56), and Ma06_g12700 from a long SOG1-like lineage found only in monocots |
| Where did SOG1 come from? | No SOG1-like gene in green algae. Moss and liverwort genes branch off **before** SOG1 split from its sister NAC group (NAC073/075/099), so they are ancestral copies of both, not strict SOG1 orthologs. The SOG1 and NAC044/085 subgroups are both present from *Amborella* onwards |
| Is the ATM/ATR switch kept? | Tails hold **146 SQ/TQ sites vs 92** expected from their composition (p < 0.001); 21 of 45 carry an SQ/TQ cluster; banana SOG1 keeps all 5 Arabidopsis SOG1 sites in the same alignment positions |
| Folded domain, flexible tail? | Median 88% of tail residues predicted disordered vs 0% of the NAC domain; matches AlphaFold confidence (ρ = −0.84 for Arabidopsis SOG1, −0.79 for banana SOG1) |
| Which part is most constrained? | Median dN/dS **0.19** in the NAC domain vs **0.40** in the tail (126 gene pairs, tail higher in 98%, p < 0.001): purifying selection on both, strongest on DNA binding |
| Which genes would MaSOG1 control? | The SOG1 site CTT(N)7AAG is enriched within 500 bp of the start of reported SOG1 target genes in Arabidopsis (**3.3×**, p < 0.001, positive control) and of their banana orthologs (**3.0×**, p < 0.001); a control motif with the same letters is not enriched |

## Pipeline

| Step | What | Tools |
|---|---|---|
| 1 | Download proteomes and CDS of 19 plant genomes (17 available), Pfam NAM domain; MD5 fingerprints | Ensembl Plants, InterPro |
| 2 | NAC family in every genome; SOG1-clade candidates by Arabidopsis matches | pyhmmer (hmmsearch, phmmer), Pfam gathering thresholds |
| 3 | Gene structure, Ensembl Compara cross-check, promoters, AlphaFold models | Ensembl REST, UniProt, AlphaFold DB |
| 4 | Alignment, ML tree, final clade call from the tree | MAFFT L-INS-i, trimAl, IQ-TREE (ModelFinder: Q.PLANT+I+G4, 1000 UFBoot) |
| 5 | SQ/TQ sites vs 10,000 shuffles, SQ/TQ cluster domains, disorder | numpy, metapredict |
| 6 | Pairwise dN/dS, NAC domain vs tail; Ks of banana paralogs | Biopython (Nei–Gojobori) |
| 7 | SOG1 motif in target vs 1,000 random promoters, with a composition control | permutation and Fisher tests |
| 8 | Interactive report | `docs/` |

```bash
pip install -r requirements.txt   # plus mafft, trimal and IQ-TREE on the PATH
python run_all.py                 # full run (needs internet, ~15 min)
python run_all.py --offline       # steps 4-8 from the data in data/derived/
python -m pytest -q
```

The GitHub workflow (`.github/workflows/pipeline.yml`) runs everything end to end and saves the results here.

## Checks built in

- NAC family sizes vs published genome-wide studies: Arabidopsis 114 (105 published), banana 181 (167 published); rice differs (100 vs 151) because the published count used a different gene annotation
- Arabidopsis SOG1, NAC044 and NAC085 must be recovered, and nothing else
- *Chlamydomonas* as a negative control (no NAC genes at all)
- SOG1 orthologs compared with Ensembl Compara: 16 of 18 found; the 2 missed are incomplete gene models (rice SOG1 is a 246-residue fragment in Ensembl)
- Every coding sequence must translate exactly to its protein before dN/dS (45 of 45 did)
- Promoter test run first on Arabidopsis, where SOG1 targets are known from ChIP-seq; the 500 bp window was set before the run
- 18 tests: unit tests on small inputs with known answers, plus checks on the output

## Limits

Sequence and motif analyses predict function; they do not show it. A motif in a promoter is not binding. Copy numbers depend on each
genome's gene models. Whether MaSOG1 rescues an Arabidopsis *sog1* mutant is a lab experiment (phase 3, listed in the report).

## Data

Ensembl Plants, InterPro/Pfam, UniProt and AlphaFold DB. `data/derived/` holds only small extracts of these public resources.

Mankirat Singh · [Portfolio](https://bymankirat.github.io/myportfolio/) · [GitHub](https://github.com/bymankirat)
