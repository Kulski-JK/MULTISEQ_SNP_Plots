# MULTISEQ_SNP_Plots

This R script performs multi-sequence SNP density analysis across genomic regions and generates pairwise SNP density plots, text summaries, and Excel tables. It is optimized for HLA haplotypes or other genomic regions of interest.

## Folder Structure

- `data/` — Input SNP tab-delimited files generated from pairwise MAUVE alignments (Darling et al 2010).
- `scripts/` — R script `multiseq_snpplots.R` for multi-sequence SNP analysis.
- `output_examples/` — Example outputs including PNG plots, TXT/SNP bin tables, and Excel files with clickable summaries.
- `.gitignore` — Prevents temporary, backup, and output files from being tracked by Git (optional).

## Features

- Reads SNP data from tab-delimited input files generated from pairwise MAUVE alignments.
- Automatically detects reference and target sequences.
- Computes SNP density per kilobase and generates histograms.
- Saves output as pairwise SNP density plots (`.png`), SNP bins per kilobase (`.txt`, `.xlsx`), and a summary Excel workbook with clickable links.
- Organizes results into an output folder automatically (`*_Rplots`).
- Designed to handle multiple pairwise comparisons in one run.

## Instructions

1. Open `scripts/multiseq_snpplots.R` in R or RStudio.
2. Ensure the required packages are installed: `ggplot2`, `openxlsx`, `dplyr`.
3. Adjust the file path in the script if necessary.
4. Run the script.
5. Outputs will be saved in a folder automatically named after the input file with `_Rplots` appended (e.g., `1v21.23.45.48.68-SNPs_Rplots`).

## Dependencies

- R (≥ 4.0)
- RStudio (≥ 2024.12.0+467)
- Packages: `ggplot2`, `openxlsx`, `dplyr`

## License

This project is licensed under the MIT License — see the [LICENSE](./LICENSE) file for details.

## References

Darling, A.C.E., Mau, B., Blattner, F.R., Perna, N.T. (2010). Mauve: multiple alignment of conserved genomic sequence with rearrangements. *Genome Research*, 14(7), 1394–1403.