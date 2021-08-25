# pharmvip-hla

Workflow of HLA module.
## Setup

### Dependencies
*   Python 3.8+.
*   [pandas 1.2.4](https://pandas.pydata.org/)

## Usage 

Create merge and summary file of 3 calculate diplotype HLA software (ATHLATES,HLAHD,KOURAMI)

```shell
python HLA_final2.py -i result_ATHLATES.txt,result_HLAHD.txt,result_KOURAMI.txt
```
Prepare result of HLA file for pharmvip-guideline.

```shell
python diplotype_HLA.py HLAMergefile.txt path/to/hla_resource sampleId
```
