# High-dimensional liver radiomics features

This repository is currently under construction and will be updated alongside our paper. Please check our paper: 

[Genetic architecture of high-dimensional liver radiomic phenotypes and their role in common metabolic diseases](https://www.medrxiv.org/content/10.64898/2026.05.19.26353617v1) 


## Overview
<img width="1037" height="685" alt="workflow" src="https://github.com/user-attachments/assets/fba33ce5-6f01-4586-87ca-581b983ca748" />




## Data
This study includes 200+ DL-derived liver radiomics features across 3D shape, first-order, and texture metrics for ~40,000 UK Biobank participants.

All liver MRI GWAS summary data will be publicly available at an appropriate time.

Main results are provided in the paper and its supplementary materials. We will add the necessary info to guide the MRI features of importance and the necessary data (eg the phenotypical correlation between QC-ed MRI features) that is anticipated to be useful in certain analyses.


## Code
Custom analysis scripts are now available at [code](code/). The analysis structures are summarized at [Code Overview](Code_overview.txt).

Details of the software and hardware dependencies/versions are provided in [System Requirements](System_requirements.txt).

## How to run the analysis
1. **Identify the analysis.** Decide which analysis from our [manuscript](https://www.medrxiv.org/content/10.64898/2026.05.19.26353617v1) you wish to run.
2. **Locate the script.** Use [Code overview](Code_overview.txt) to find the corresponding script in the [code](code/) folder.
3. **Check dependencies.** Software and package versions are listed in [System requirements](System_requirements.txt).
4. **Run the script.** Each script is self-contained, with detailed guidance provided as in-script comments. Update the file paths to point to your input data. 

## Citation

If you use the code in this repository, or build on the analysis design or pipelines described in our paper, please cite:

Tian H, Kamineni M, Truong B, et al. Genetic architecture of high-dimensional liver radiomic phenotypes and their role in common metabolic diseases. medRxiv. 2026.

```bibtex
@article{tian2026MRI,
  title   = {Genetic architecture of high-dimensional liver radiomic phenotypes and their role in common metabolic diseases},
  author  = {Tian, Haodong and Kamineni, Manish and Truong, Buu and others},
  journal = {medRxiv},
  year    = {2026}
}
```

## Questions
Please open an issue for any question.
