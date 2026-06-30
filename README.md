# TenFuse

TenFuse is a MATLAB implementation of the blind hyperspectral and multispectral image fusion framework proposed in:

> Y. Gao, M. K. Ng, and C. Cui, "Blind Hyperspectral and Multispectral Images Fusion: A Unified Tensor Fusion Framework from Coupled Inverse Problem Perspective."

The code simulates low-spatial-resolution hyperspectral images (LR-HSI) and high-spatial-resolution multispectral images (HR-MSI) from the Washington DC Mall hyperspectral image, and then reconstructs the target high-spatial-resolution hyperspectral image (HR-HSI) using the proposed tensor fusion model and ADMM-based optimization algorithm.

## Repository Structure

```text
TenFuse/
|-- DCmall_synthetic_fusion.m      # Main demo script
|-- DC_mall.mat                    # Washington DC Mall hyperspectral data
|-- SRF_matrix.mat                 # Spectral response matrix
|-- Initialization/                # Initialization methods
|-- My_Solver/                     # TenFuse optimization solver
|-- functions/                     # Utility functions for tensor and image operations
|-- quality_indices/               # Quantitative evaluation metrics
```

## Requirements

This project is implemented in MATLAB. The code uses standard MATLAB functions and functions from the Image Processing Toolbox, such as `fspecial` and `psf2otf`.

Recommended environment:

- MATLAB R2020a or later
- Image Processing Toolbox

## Quick Start

1. Clone this repository:

```bash
git clone https://github.com/gaoYn-start/TenFuse.git
cd TenFuse
```

2. Open MATLAB and set the repository folder as the current working directory.

3. Run the main script:

```matlab
DCmall_synthetic_fusion
```

The script will:

- Load the Washington DC Mall hyperspectral image and spectral response matrix;
- Generate simulated LR-HSI and HR-MSI observations;
- Initialize the spatial and spectral degradation matrices;
- Run the TenFuse ADMM solver;
- Produce the reconstructed HR-HSI result `SRI_tenfuse` and record the running time `time_tenfuse`.

## Data

The file `DC_mall.mat` contains the Washington DC Mall hyperspectral image used as the reference image in the demo. The file `SRF_matrix.mat` stores the spectral response matrix used to generate the simulated multispectral observation.

Note that `DC_mall.mat` is a relatively large file. GitHub may display a large-file warning when downloading or cloning the repository.

## Main Output

After running `DCmall_synthetic_fusion.m`, the main output variables are:

- `SRI_tenfuse`: reconstructed high-spatial-resolution hyperspectral image;
- `time_tenfuse`: computational time of the TenFuse algorithm.

## Citation

If you use this code in your research, please cite the corresponding paper:

```bibtex
@article{gao_tenfuse,
  title   = {Blind Hyperspectral and Multispectral Images Fusion: A Unified Tensor Fusion Framework from Coupled Inverse Problem Perspective},
  author  = {Gao, Ying and Ng, Michael K. and Cui, Can},
  journal = {To appear},
  year    = {2026}
}
```

The source code for reproducing our experiments is publicly available at:

```text
https://github.com/gaoYn-start/TenFuse
```

## Contact

For questions about the code, please contact the repository maintainer.
