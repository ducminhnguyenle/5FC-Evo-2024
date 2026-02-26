# 5FC-Evo-2024

## Overview

This repository accompanies our paper:

Phan-Canh T*, Nguyen-Le D-M*, Luu P-L, Khunweeraphong N, Kuchler K. Rapid *in vitro* evolution of flucytosine resistance in *Candida auris*. mSphere. 2025;10(4):e0097724. doi:10.1128/msphere.00977-24

> <https://journals.asm.org/doi/10.1128/msphere.00977-24>

## 📝 Analysis Workflow

- Main workflow, data analysis and visualization source code availability: <https://github.com/ducminhnguyenle/5FC-Evo-2024/>
- SNP variant calling workflow using bulk RNAseq data: [SNP calling workflow](./SNP_calling/)
- SNP variant annotation and profiling from RNAseq datasets: [SNP annotation](./SNP_annotation/)
- WGS data analysis and variant annotation: [WGS analysis](./WGS)
  - Cut-off **80%** ALT reads: [WGS analysis](./WGS/alt_AD_80/)
  - Cut-off **40%** ALT reads: [WGS analysis](./WGS/alt_AD_40/)
- CNV calling from WGS data: [CNV analysis](./WGS/CNV/)
- Docking analysis for single and multiple ligands simultaneously: [Docking](./docking_analysis/)

## Citation

If you find this repository helpful, please cite our paper:

```bibtex
@article{PhanCanh2025,
  author  = {Phan-Canh, Trinh and Nguyen-Le, Duc-Minh and Luu, Phuc-Loi and Khunweeraphong, Narakorn and Kuchler, Karl},
  title   = {Rapid in vitro evolution of flucytosine resistance in Candida auris},
  journal = {mSphere},
  year    = {2025},
  volume  = {10},
  number  = {4},
  pages   = {e00977-24},
  doi     = {10.1128/msphere.00977-24},
  url     = {https://journals.asm.org/doi/10.1128/msphere.00977-24},
  note    = {Trinh Phan-Canh and Duc-Minh Nguyen-Le contributed equally to this work}
}
```

## License

This repository is licensed under the [MIT License](LICENSE).
