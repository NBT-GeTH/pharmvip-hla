# pharmvip-hla

Workflow of HLA module.
## Setup

### Dependencies
*   Python 3.8+.
*   [pandas 1.2.4](https://pandas.pydata.org/)

## Usage 

Create output files from result summarization of three HLA genotyping software (ATHLATES,HLA-HD,KOURAMI)

```shell
python HLA_final2.py -i result_ATHLATES.txt,result_HLAHD.txt,result_KOURAMI.txt
```
Prepare HLA result file for pharmvip-guideline module.

```shell
python diplotype_HLA.py HLAMergefile.txt path/to/hla_resource sampleId
```
