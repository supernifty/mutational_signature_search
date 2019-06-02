## Mutational Signature Search

Calculate mutational signatures over many samples with many experimental parameters.

## Installation
On spartan:
```
module load Python/3.6.4-intel-2017.u2
module load cURL/7.60.0-spartan_gcc-6.2.0
```

```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage

* plot_hist: plot evolution of signatures across a parameter for a single sample

## Inputs
* VCF files
* AF points
* DP points
* filename template
