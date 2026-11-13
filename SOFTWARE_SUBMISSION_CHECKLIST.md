# Software documentation checklist

Documentation follows the [Nature code-publication guidance](https://media.nature.com/full/nature-cms/documents/GuidelinesCodePublication.pdf): installation, tested dependencies, inputs, demo commands, expected output, runtime and access conditions. The target journal's requirements take precedence.

| Item | Location/status |
|---|---|
| Code purposes and manuscript mapping | CODE_MAP.md/CSV |
| Field definitions and units | DATA_DICTIONARY.md; DATA_SCHEMA.csv |
| Pinned tested versions | requirements.txt; Python 3.12 |
| Installation and demo | README.md; `python run_all.py` |
| Expected results | fPC1=81.1%, fPC2=14.9%, silhouette 0.790; 156 RV pairs and 720 permutations/pair |
| Numerical verification | VALIDATION.md and generated verification reports |
| Runtime | Per-module times in generated results/run_times.json; no GPU required |
| Confidential verification | Separate reviewer dataset and explicit command |
| Version | Local release candidate 1.0.0; cite final published version/commit |
| Identity | Anonymous code/data filenames and no personal file paths or individual public records |
| Licensing | LICENSE_STATUS.md; rights-holder selection pending |
| Availability | Actual repository URL will be entered after upload |

The optional local figure-audit hook is not required to run the public code. No original instrument or unrelated exploratory workflow is represented as a manuscript analysis.
