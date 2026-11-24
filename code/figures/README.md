# Figure builders

Run the builders from the code repository root after staging the preserved outputs:

```sh
python scripts/stage_source_data.py --source "../source data" --code-root .
```

| Folder | Contents |
| --- | --- |
| `main/` | Builders for Figures 1-5 |
| `supplementary/` | Builder for Supplementary Figure S3 |
| `graphical_abstract/` | Python builder for the graphical abstract |
| `shared/` | R plotting functions, colours, number labels and data access |
| `data/` | Supplied numerical manifest |

For example:

```sh
Rscript figures/main/make_fig01_construction.R
Rscript figures/supplementary/make_figS_reuse_native_agreement.R
python figures/graphical_abstract/build_graphical_abstract.py
```

These builders write figures beside their scripts. Supplementary Figures S1 and S2 are produced by `analysis/sensitivity/report_patient_omission.py`, which also reports the patient-omission analyses.

The R builders require ggplot2, patchwork, jsonlite, scales and ragg. The graphical abstract requires ReportLab, Pillow, Poppler and Arial font files at the path specified in its builder. Its diagrams are drawn directly from vector shapes in code.
