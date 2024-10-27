## Mutational Signature Search

Calculate mutational signatures over many samples with many experimental parameters.

## Installation
```
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Usage

### Step 1: Generate a table of signatures at different settings
```
python mutational_signature_search/search.py --signatures sig_def --genome genome.fa --vcfs input_vcfs --dps depths --afs afs --caller strelka --use_bam_depth --tags wes-strelka-pass > search.out
```

*Examples*

TCGA
```
python mutational_signature_search/search.py --genome ../somatic_pipeline_hg38/pipeline/assets/Homo_sapiens_assembly38.fasta --signatures ../src/mutational_signature/data/signatures_cosmic_v3.4_sbs.hg38.txt --vcfs example/maf.gz --maf_sample TCGA-AA-3966-01A-01D-1981-10 --caller maf --dp 10 50 100 --af 0.1 0.3 0.5 --use_bam_depth t_depth --use_alt_depth t_alt_count
```


Some notes:
* sig_def can be taken from the mutational_signatures package
* depths is a list of depths to test e.g. 50 100 150 200
* same for afs
* use_bam_depth is a flag to use if GATK's bamdepth tool has been run on the vcf

The tool can be slow to run, so consider running the tool multiple times with different combinations of depths and afs, then combine the files together.

This is the main computationally expensive step, once this is complete, there are different ways to visualise the data.

Expectation for plotting is a tsv with the following columns:
* DP
* AF
* Sample
* Tags
* Error
* Variants
* Signature values

### plot evolution of signatures across a parameter for a single sample
Example command line:
```
python software/mutational_signature_search/mutational_signature_search/plot_heat.py --data plot_search --x DP --y AF --highlight SBS36 --target out/sample.heat.png --samples sample_name
```

### plot individual signature across multiple parameters for a single sample (heat map)

## Inputs
* VCF files
* AF points
* DP points
* filename template
