# EukAsmChecker

As more eukaryotic genomes are sequenced and assembled, ensuring assembly quality becomes increasingly important. Contamination introduced during sampling, cultivation, or sequencing, particularly from prokaryotes, may escape routine quality control. Along with assembly errors, these contaminants can compromise downstream gene prediction, functional annotation, and comparative genomics.

We developed EukAsmChecker to help researchers inspect and filter genome assemblies through an interactive interface. It combines sequence length, reads coverage, and classification labels from [EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner) (ECM), a tool that uses language models to classify sequences with high accuracy. Users can explore sequence distributions, adjust filtering thresholds, choose which labels to retain, and export filtered FASTA files and summary tables.

Both web and command line versions are available to make EukAsmChecker accessible and convenient to use across different workflows. README is also available in [Chinese](README.zh-CN.md).

## Web version

[Open the web app](https://westraingroup.github.io/EukAsmChecker/) to review assemblies without installing software.

1. Select uncompressed FASTA, coverage TSV and ECM TSV files. Multiple samples are supported.
2. Check the sample-to-file assignments and select the column names.
3. Review the plots and filters, then confirm each sample.
4. Click **Save results** to download one ZIP containing filtered FASTA, sequence tables, a project summary and decisions.json. Examples contain sequence statistics and export tables and decisions.

Refreshing the web page starts a new session. Each visit to Example starts from the defaults. Save results before refreshing or closing the page. To resume your own saved review, select the same input files and use **Import review (.json)** to load decisions.json from the result archive. Column settings apply to every sample.

## Command-line version

### Install

Python 3.8 or newer and the cp command are required on Linux or macOS.

Install with Conda (Linux or macOS):

```bash
conda create -n eukasmchecker --override-channels -c westraingroup -c conda-forge eukasmchecker
conda activate eukasmchecker
eukasmchecker --help
```

Alternatively, install from GitHub:

```bash
python -m pip install git+https://github.com/WeStrainGroup/EukAsmChecker.git
```

### Prepare

```bash
eukasmchecker prepare \
  --sample sample_A \
  --fasta assembly.fasta \
  --ecm classification.tsv \
  --coverage coverage.tsv \
  --seq-id scaffold_id --seq-cov mean_depth \
  --out project/sample_A/sample_A.eac.json
```

FASTA IDs use the first whitespace-delimited field after >. All three inputs must contain exactly the same IDs, once each. Coverage values must be finite and nonnegative. Lengths in the ECM table must match FASTA.

The coverage table has a header. Select the ID and coverage columns with --seq-id and --seq-cov.

```tsv
scaffold_id	mean_depth
contig_1	88.8368
contig_2	12.5
contig_3	7.25
contig_4	0
```

The CLI reads contig_id, length_bp and label_euk from the ECM table. The web app allows you to select equivalent columns. The label euk is normalized to eukaryota. Additional columns are accepted.

```tsv
contig_id	length_bp	label_euk
contig_1	25000	eukaryota
contig_2	1500	other
contig_3	750	eukaryota
contig_4	499	NA
```

### Review

Try the review interface with the [interactive example](https://westraingroup.github.io/EukAsmChecker/?demo=1).

Click a point to toggle its selection. Drag a rectangle to select its sequences, or to deselect them when the whole region is already selected. Other selections are preserved. Use **Remove selected** or Delete to apply the exclusions. **Restore manual removals** restores those sequences to threshold-based filtering. Manual removals are saved in decisions.json and applied to FASTA and tables.

```bash
eukasmchecker review project/ --threads 8
# Create a portable offline review page
eukasmchecker review project/ --out review.zip --threads 8
```

A directory is searched recursively for .eac.json files. Direct review opens a local browser session and saves progress to .eukasmchecker/<dataset_id>/decisions.json. Keep the command running while reviewing.

For offline review, download and extract the complete ZIP, then open index.html. Keep its data folder beside it. Save decisions downloads decisions.json.

### Apply

```bash
eukasmchecker apply decisions.json --out project/ --threads 8
```

Apply verifies the saved review, source files and retained counts. Results can be written beside the prepared inputs or into a new output directory. Existing result files are protected from overwrite.

```text
project/
└── sample_A/
    ├── sample_A.eac.json
    ├── sample_A.fasta
    └── sample_A.sequences.tsv
```

FASTA preserves the original retained records byte for byte, including headers, wrapping, order and line endings. The TSV contains every original sequence. For the four input rows above, default thresholds and all labels selected produce the following table.

```tsv
seq_id	length_bp	coverage	ecm_label	filter_status	filter_reason
contig_1	25000	88.8368	eukaryota	retained	pass
contig_2	1500	12.5	other	retained	pass
contig_3	750	7.25	eukaryota	removed	low_coverage
contig_4	499	0	NA	removed	short_length;low_coverage
```

Removal reasons are short_length, low_coverage, high_coverage, excluded_label and manual_exclusion. Multiple reasons are joined with semicolons. Excluding a sample produces an empty FASTA and excluded_sample for every sequence.

## Filtering rules

The web and command-line versions use the same rules.

Scatter points use a reproducible random drawing order across labels (seed 6311).

| Setting | Default and behavior |
|---|---|
| Minimum length | 500 bp. |
| Minimum coverage | `round(max(reference median / 25, 10))`. The reference set is original sequences labelled `eukaryota` with length >10 kb. Set a minimum manually if this set is empty. |
| Maximum coverage | The original maximum, rounded upward; no sequences are removed by default. |
| Included labels | All labels are selected initially; any nonempty subset can be retained. |
| Manual removals | None initially. Selected sequence IDs remain excluded until restored. Reset defaults also clears manual removals. |
| Retain sample | On by default. Turning it off removes every sequence in that sample. |

Retention uses the length minimum, both inclusive coverage bounds, selected labels and saved manual removals. Coverage is the per-sequence mean from the input table. Reference sequences receive equal weight, including zero coverage. Thresholds are integers, with half values rounded upward; filtering and output tables preserve original coverage precision. For example, 84.9343× is below an 85× minimum even when displayed as 85×.

[EukContigMiner](https://github.com/WeStrainGroup/EukContigMiner#metrics) (ECM) is less reliable for short contigs. For 500–1,000 bp contigs, the test F1 is approximately 0.9 at 1% eukaryotic prevalence. Some short contigs within the main assembly cluster may therefore be labelled `other` through classification error; assess their coverage and cluster position together with the label. ECM does not classify sequences below 500 bp, which retain `NA`.

## Files and saved reviews

Prepared CLI inputs and decisions record absolute compute-platform paths. Keep the prepared files and original FASTA at those paths. After offline review, return decisions.json to that platform before applying it. The original coverage and ECM tables are no longer required after preparation.

Browser decisions refer to browser-prepared files and can be restored after selecting the same files in the web app. CLI apply uses its own prepared-file references. Browser project capacity depends on available memory and storage.

## License

Copyright © 2026 WeStrainGroup.

**CC BY-NC-SA 4.0**. See [LICENSE](LICENSE). The bundled js-sha256 dependency retains its [MIT license](web/src/vendor/LICENSE.js-sha256.txt).
