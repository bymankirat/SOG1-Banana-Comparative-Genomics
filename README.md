# SOG1 in banana: phase 2 of my *Ma-SOG1* project

**Live report:** https://bymankirat.github.io/SOG1-Banana-Comparative-Genomics/

SOG1 (SUPPRESSOR OF GAMMA RESPONSE 1) is the plant counterpart of p53. After DNA damage the ATM and ATR kinases phosphorylate
SQ/TQ sites in its tail, and SOG1 then switches on DNA repair, stops the cell cycle or triggers cell death.

In 2023, during research training at the National Agri-Food Biotechnology Institute (NABI), Mohali, I analysed the banana gene
*Ma-SOG1* in silico and cloned it into a plant expression vector for overexpression in *Arabidopsis*. That lab work belongs
to the NABI group and is not included here. This repository is the computational next step, built only on public genomes.

## Questions

1. **How many SOG1 genes does banana have?** Banana went through three whole-genome duplications, so one gene in Arabidopsis may be several in banana.
2. **Where do they sit in plant evolution?** A maximum-likelihood tree of the whole SOG1 clade, from moss to banana.
3. **Is the ATM/ATR switch conserved?** SQ/TQ sites in the tail, tested against chance; disorder vs AlphaFold.
4. **Which parts are most constrained?** dN/dS of the DNA-binding domain vs the tail.
5. **Which banana genes would MaSOG1 control?** The SOG1 binding motif CTT(N)7AAG in promoters of known targets: Arabidopsis as positive control, then banana orthologs.

## Pipeline

| Step | What | Tools |
|---|---|---|
| 1 | Download proteomes and CDS of 19 plant genomes, Pfam NAM domain (MD5 fingerprints) | Ensembl Plants, InterPro |
| 2 | NAC family in every genome; SOG1 clade by best Arabidopsis match | pyhmmer (hmmsearch, phmmer), Pfam gathering thresholds |
| 3 | Gene structure, Ensembl Compara cross-check, promoters, AlphaFold models | Ensembl REST, UniProt, AlphaFold DB |
| 4 | Alignment and ML tree | MAFFT L-INS-i, trimAl, IQ-TREE 2 (ModelFinder, 1000 UFBoot) |
| 5 | SQ/TQ sites vs 10,000 shuffles, SQ/TQ cluster domains, disorder | numpy, metapredict |
| 6 | Pairwise dN/dS, NAC domain vs tail; Ks of banana paralogs | Biopython (Nei–Gojobori) |
| 7 | SOG1 motif in target vs random promoters, with a composition control | permutation and Fisher tests |
| 8 | Interactive report | `docs/` |

```bash
pip install -r requirements.txt   # plus mafft, trimal, iqtree2 on the PATH
python run_all.py                 # full run (needs internet)
python run_all.py --offline       # steps 4-8 from the data in data/derived/
python -m pytest -q
```

The GitHub workflow (`.github/workflows/pipeline.yml`) runs everything end to end and saves the results here.

## Checks built in

- NAC family sizes compared with published genome-wide studies (Arabidopsis, rice, banana)
- Arabidopsis SOG1, NAC044 and NAC085 must be recovered, and nothing else
- Green alga *Chlamydomonas* as a negative control (SOG1 is a land-plant gene)
- Clade membership compared with Ensembl Compara gene trees
- Every coding sequence must translate exactly to its protein before it is used for dN/dS
- Promoter test run first on Arabidopsis, where SOG1 targets are known from ChIP-seq
- Unit tests on small inputs with known answers

## Limits

Sequence and motif analyses predict function; they do not show it. Whether MaSOG1 rescues an Arabidopsis *sog1* mutant is a lab experiment (phase 3, listed in the report).

## Data

Ensembl Plants, InterPro/Pfam, UniProt and AlphaFold DB. `data/derived/` holds only small extracts of these public resources.

Mankirat Singh · [Portfolio](https://bymankirat.github.io/myportfolio/) · [GitHub](https://github.com/bymankirat)
