# Carbon-A Gene Prediction + BUSCO (Colab)

Batch *ab initio* gene prediction on Google Colab using
[Carbon-A-1.2B](https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B) — Hugging Face Biology's
1.2-billion-parameter DNA language model that predicts, for every base on both strands, the
probability that it is protein-coding — with a gene-grammar Viterbi decoder that turns those
probabilities into gene models, bigWig probability tracks, and built-in **BUSCO** completeness
assessment. Point it at a genome (or a whole folder of them) and it returns per-genome GFF3 + CDS +
protein FASTA, coding-probability tracks, a prediction stats table, and BUSCO scores, all from a single
notebook with form-field inputs.

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/amiralito/CarbonA_colab/blob/main/CarbonA_gene_prediction.ipynb)

## Pipeline

```
genome FASTA  →   Carbon-A-1.2B (GPU)                 →   gene-grammar Viterbi (CPU)       →   GFF3 + CDS +   →   BUSCO
(one or many)     98,304-bp windows, 6-mer tokens         ATG … first in-frame stop,            protein FASTA      (proteins and/or
                  per-base CDS probability, both          GT-AG / GC-AG / AT-AC introns,        + bigWig tracks    genome mode)
                  strands; overlaps averaged;             frame kept across junctions
                  checkpointed → resumable
```

## About the gene models — read this first

The public release of Carbon-A is the **annotator only**: it outputs per-base CDS probabilities, not
genes. The [technical report](https://huggingface.co/buckets/HuggingFaceBio/Carbon-A-Technical-Report)
turns them into gene models with three further stages — a small learned *transition network*, a
gene-grammar *Viterbi* decoder, and a learned *confidence network* — which are **not released**.

Cell 8 re-implements the Viterbi stage as the report describes it:

- a gene begins at an **ATG** and ends at the **first in-frame stop** under the selected genetic code
  (stops spanning a splice junction included); the start codon is never split by an intron;
- an intron runs from a donor to a compatible acceptor — **GT…AG, GC…AG or AT…AC** — with fixed
  log-priors from their genomic frequencies (0.9935 / 0.0056 / 0.0009), and the reading frame is
  carried across the junction;
- emissions are **log p** for CDS bases and **log(1 − p)** for intergenic and intronic bases;
- records may start or end inside a gene — such models are kept and flagged `partial=`.

What differs from the report: the learned transition network is replaced by two fixed, tunable
penalties (`GENE_PENALTY`, `INTRON_PENALTY`) plus a minimum intron length, and the trained confidence
score is replaced by a simple **gene score** (mean CDS probability × mean intronic non-coding
probability). Like the report: **one isoform per locus, CDS only (no UTRs)**.

On a 120-kb stretch of *Arabidopsis thaliana* chromosome 1 (31 RefSeq loci), the raw model reached a
nucleotide-level CDS F1 of 0.976 and the decoder reproduced the exact RefSeq CDS structure of one of
the isoforms at 26 of 31 loci; the misses were one locus split into two genes and four loci with a missed, extra or mis-bounded
exon (two of them via spurious 3–4-bp micro-exons). Results were stable for penalties between 0 and 8.
Treat this as a sanity check, not a benchmark.

## What it does

- Runs Carbon-A over **one or many genomes** (recursive folder scan, upload, or URL); FASTA may be
  gzipped, soft-masked genomes need no preprocessing.
- Loads the model straight into **Colab's own PyTorch** (`transformers ≥ 4.56`, custom model code from
  the commit-pinned snapshot) — no separate environment, no restart.
- Downloads the **weights once** (4.6 GB) into `WORKDIR/carbonA_model/<commit>/` and resolves
  `MODEL_REVISION` to a commit hash recorded in `run_config.json`.
- **Verifies** the setup on your GPU: the built-in vectorised 6-mer encoder against the model's own
  tokenizer, finite output on a full window of your genome, fp16 vs fp32 agreement on a T4, and a
  **measured run-time estimate** for the whole batch.
- Uses the report's genome-scale windowing: 98,304-bp windows overlapping by 6,144 bp, last window
  right-aligned to the sequence end, predictions averaged where windows overlap (`49152` reproduces
  the model card's 50 % overlap).
- **Checkpoints** window predictions (int8 logits, ~1 byte per bp) every few minutes, so long genomes
  survive disconnects: re-open the run folder with `RESUME_RUN` and only missing windows are predicted.
- **Decodes** gene models on the CPU from the saved windows — re-run cell 8 alone to try other
  penalties, genetic code or filters without touching the GPU.
- Writes **GFF3** (gene, mRNA, exon, CDS with phases; score + `partial=` attributes), **CDS + protein
  FASTA**, and **bigWig** CDS-probability tracks per strand for IGV / JBrowse.
- Supports alternative **nuclear genetic codes** (1, 6, 10, 12, 15, 26, 29, 30) — e.g. ciliates.
- **Gzips** outputs into `gff/`, `cds/`, `protein/`, `tracks/`, `log/`, `busco/` inside a
  **timestamped run folder** so consecutive runs never overwrite each other.
- Builds a per-genome **stats table** and scores completeness with **BUSCO** in proteins and/or genome
  mode, compiling a summary table that accumulates across lineages.
- All code cells are collapsed to **form fields**; mount Google Drive to persist inputs, weights and
  outputs.

## Requirements

- Google Colab with a **GPU runtime** (`Runtime → Change runtime type → GPU`).
- **A100 / L4 / H100 (bf16) strongly recommended.** A free **T4** works in fp16 but is several times
  slower; cell 6 prints the measured time per window and an estimate for your genomes. A CPU runtime
  is only usable for tiny tests (~10 min per 100 kb).
- One full window needs ~6 GB of GPU memory at batch 1.
- ~4.6 GB of Drive space for the cached weights, plus ~1 byte per bp of genome for the saved windows.
- Host RAM for decoding: ~17 bytes per bp of the **largest chromosome** (use a high-RAM runtime for
  chromosomes beyond ~500 Mb).
- Genome FASTA only (`.fa/.fasta/.fna/.fas`, optionally `.gz`). No RNA-seq or protein evidence.

## Quick start

1. Open the notebook in Colab and select a **GPU runtime**.
2. Run **cell 1** (runtime check) and **cell 2** (install).
3. **Cell 3** — optionally mount Drive, set `WORKDIR`, and stamp this run's folder (or resume one).
4. **Cell 4** — set prediction options and choose the genome input mode.
5. **Cell 5** — collect the genome(s) and scan their sequences.
6. **Cell 6** — download and verify the model; read the run-time estimate.
7. **Cell 7** — run Carbon-A (GPU); **cell 8** — decode gene models; **cell 9** — view the stats table.
8. **Cells 10–13** — configure and run BUSCO, then compress the results.

## Notebook cells

| Cell | Step |
|------|------|
| 1 | Check runtime — GPU, memory, compute capability |
| 2 | Install dependencies (`transformers` if < 4.56, `pyBigWig`) and write `carbon_core.py` |
| 3 | Storage — mount Drive & stamp (or resume) a run folder |
| 4 | Prediction parameters & genome input |
| 5 | Collect genome(s) and scan sequence lengths |
| 6 | Download Carbon-A (commit-pinned) and verify it on this GPU; run-time estimate |
| 7 | Run Carbon-A on all genomes → checkpointed window predictions in `windows/` |
| 8 | Decode gene models → gzipped `gff/ cds/ protein/` + bigWig `tracks/` |
| 9 | Per-genome prediction stats table → `prediction_summary.tsv` |
| 10 | BUSCO options (mode toggles + lineage) |
| 11 | Install BUSCO (+ genome-mode tools if enabled) |
| 12 | Run BUSCO + compile `busco_summary.tsv` |
| 13 | Compress & save BUSCO run directories into `carbonA_predictions/busco/` |

## Inputs & parameters

**Model (cell 4)**
- `MODEL_REVISION` — branch or commit of `HuggingFaceBio/Carbon-A-1.2B`; resolved to a commit hash in
  cell 6. Paste that hash to repeat a run exactly.
- `PRECISION` — `auto` = bf16 on Ampere and newer, fp16 on T4 / V100, fp32 on CPU.
- `WINDOW_OVERLAP` — `6144` (report, default) or `49152` (model card's `infer_genbank.py`, ~2× GPU time).
- `BATCH_SIZE` — windows per forward pass; `auto` = 1 / 2 / 4 for < 40 / 40 / 80 GB. Halved
  automatically on OOM.
- `MIN_RECORD_LENGTH` — skip sequences shorter than this (default 1,000 bp; the report's pipeline uses
  100,000 bp, and predictions on very short contigs are less reliable).
- `CHECKPOINT_MINUTES` — how often window predictions are flushed to the run folder.

**Decoder (cell 8)** — re-run cell 8 alone to change these.
- `PREFIX` (default `CarbonA`) — gene IDs `CarbonA_<seqid>_000001`, transcripts `.1`. Empty = genome
  file name.
- `GENE_PENALTY` (8), `INTRON_PENALTY` (4) — log-penalties for opening a gene / an intron. Higher =
  fewer, longer genes / fewer introns. Coding evidence dominates; results are stable over ~0–8.
- `MIN_INTRON_LENGTH` (20) — shortest intron the decoder may place.
- `GENETIC_CODE` — NCBI table deciding stops and translation (1 standard, 6 ciliate, 10 euplotid,
  12 alternative yeast, 15 *Blepharisma*, 26 *Pachysolen*, 29 *Mesodinium*, 30 peritrich).
- `MIN_CDS_LENGTH` (96, the report's default), `MIN_SCORE` (0 = keep all), `KEEP_PARTIAL` (on).
- `WRITE_BIGWIG` — per-strand CDS-probability tracks (3 decimals).

**Batch behaviour**
- `SKIP_DONE` — skip genomes already predicted in this run folder.
- `STOP_ON_ERROR` — off by default: a failed genome is logged, listed at the end of the cell, and the
  batch moves on.

**Genome source** (`GENOME_SOURCE`)
- `drive_dir` — run on **every** FASTA found recursively under `GENOME_DIR`.
- `drive_file` — a single genome already on Drive (`GENOME_PATH`).
- `url` — download a genome from `GENOME_URL`.
- `upload` — upload from your machine.

## Outputs

```
WORKDIR/
├── carbonA_model/<commit>/                   # cached weights (4.6 GB), shared by every run
└── runs/
    └── run_20261009-142530[_tag]/            # stamped by cell 3
        ├── carbonA_predictions/
        │   ├── gff/      <genome>_gff.gff3.gz          # gene, mRNA, exon, CDS (+ score, partial=)
        │   ├── cds/      <genome>_cds.fasta.gz
        │   ├── protein/  <genome>_protein.fasta.gz
        │   ├── tracks/   <genome>_cds_prob_plus.bw, <genome>_cds_prob_minus.bw
        │   ├── windows/  <genome>/chunk_*.npz, manifest.json, DONE.json   # raw window predictions
        │   ├── log/      <genome>.carbonA.log
        │   └── busco/    <genome>_busco_{prot,genome}_<lineage>.zip + busco_summary.tsv
        ├── prediction_summary.tsv   # per-genome: size, genes, genes/Mb, partial, protein lengths, scores
        ├── run_config.json          # model commit, windowing, precision, GPU, package versions
        └── decode_config.json       # decoder settings (a re-decode adds a timestamped copy)
```

The FASTA headers carry the gene ID, location, exon count, score and any `partial=` flag. Delete
`windows/` once you are happy with the decoding if you need the Drive space back — but cell 8 can no
longer re-decode without it.

## Notes & gotchas

- **Resuming:** set `RESUME_RUN` in cell 3 to the interrupted run folder and run cells 1–7 again.
  Window size, overlap, model commit and FASTA must be unchanged — cell 7 checks and refuses mixtures.
- **Soft-masking is fine** — sequences are upper-cased before tokenization, as the model requires
  (lowercase bases would otherwise not form 6-mers). Don't hardmask: `N` can never be coding, so a gene
  is cut at an `N`.
- **Overlapping genes** on opposite strands are kept (strands are decoded independently); on the same
  strand one gene per position and one isoform per locus.
- **Scope:** eukaryotic nuclear genomes. Bacteria and viruses are outside the training data.
- **GPU OOM at batch 1** means something else holds GPU memory — restart the runtime.
- **Edit a hidden cell** → cell menu (⋮) → *Show code*.

## Citation

If you use this notebook, please cite Carbon-A and the tools it depends on.

- **Carbon-A** — HuggingFaceBio. *Carbon-A* technical report (2026) and model card,
  https://huggingface.co/HuggingFaceBio/Carbon-A-1.2B
- **GENERanno** — the GenerTeam annotator Carbon-A is derived from
  (`GenerTeam/GENERanno-eukaryote-1.2b-cds-annotator-preview`).
- **BUSCO** — Manni M, Berkeley MR, Seppey M, Simão FA, Zdobnov EM. *BUSCO update: novel and
  streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic,
  prokaryotic, and viral genomes.* Mol Biol Evol. 2021. doi:10.1093/molbev/msab199
