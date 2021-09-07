# pharmvip-hla

Workflow of HLA module.
## Setup

### Dependencies
*   Python 3.8+.
*   [pandas 1.2.4](https://pandas.pydata.org/)

## Usage 

1. Create output files from result summarization of three HLA genotyping software (ATHLATES,HLA-HD,KOURAMI)
- input : result_ATHLATES.txt result_HLAHD.txt result_KOURAMI.txt
- output : hla_summary.txt hla_results_for_diplotype.txt
```shell
python HLA_final2.py -i result_ATHLATES.txt,result_HLAHD.txt,result_KOURAMI.txt
```
2. Prepare HLA result file for pharmvip-guideline module (parameter --diplotype_hla).
- input : hla_results_for_diplotype.txt haplotype_mapping_HLA.txt
- output : diplotype_HLA.tsv

```shell
python diplotype_HLA.py hla_results_for_diplotype.txt path/to/hla_resource/haplotype_mapping_HLA.txt NameSample
```
