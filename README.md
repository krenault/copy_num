# copy_num

Phylogenetic analysis of mammalian **gene copy-number** variation (TOGA orthology) in relation to longevity traits, plus a reusable **reciprocal BLAST** tool to validate putative gains/losses.

## What's here

| Path | Description |
|------|-------------|
| `scripts/` | TOGA assembly → PCA/UMAP → elastic net / PGLMM / phylo-ranked associations → plots |
| `data/` | Species tree, metadata, orthologous copy-number matrix |
| `results/` | Example association and enrichment tables |
| `plots/` | Example figures |
| `tools/reciprocal_blast/` | API-based reciprocal BLAST validation (standalone) |

## Quick start (analysis)

Run scripts from the **repo root** (or set `COPY_NUM_ROOT`):

```bash
# Example: phenotype assembly expects optional external AnAge-style tables under data/external/
Rscript scripts/2_TOGA_PCA_orthology.R
Rscript scripts/4_PGLMM.R
```

GMT gene-set files for fgsea are read from `GSEA_DIR` or `data/gsea/` (not bundled; download MSigDB locally).

## Quick start (validation tool)

```bash
cd tools/reciprocal_blast
pip install -r requirements.txt
python validate_gene.py --gene NKG7 --target-species "Myotis davidii" --protein
```

Details: [`tools/reciprocal_blast/README.md`](tools/reciprocal_blast/README.md).

## Notes

- Personal absolute paths have been replaced with repo-relative `ROOT` / `COPY_NUM_ROOT`.
- Large genome downloads and Python virtualenvs are **not** included.
- Some early assembly scripts still expect optional external databases under `data/external/` — see script headers.

## License

MIT
