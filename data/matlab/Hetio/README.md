### Data matrices for Hetionet: https://github.com/dhimmel/hetionet.

### The matlab matrices can be imported from python.
```py
from scipy.io import savemat, loadmat
x=loadmat('hetio_AdG.mat')
```
### Matrices information
#### Adjacency matrix(i,j)=1 if ith row and jth column is positive  
```
Name      : Edge type                      : Shape                : #positive   : density
hetio_AdG : Anatomy-downregulates-Gene     : 402 x 20945          : 102240      : 0.0121
hetio_AeG : Anatomy-expresses-Gene         : 402 x 20945          : 526407      : 0.0625
hetio_AuG : Anatomy-upregulates-Gene       : 402 x 20945          : 97848       : 0.0116
hetio_CbG : Compound-binds-Gene            : 1552 x 20945         : 11571       : 3.6e-4
hetio_CcSE: Compound-causes-Side Effect    : 1552 x 5734          : 138944      : 0.0156
hetio_CdG : Compound-downregulates-Gene    : 1552 x 20945         : 21102       : 6.5e-4
hetio_CpD : Compound-palliates-Disease     : 1552 x 137           : 390         : 0.0018
hetio_CtD : Compound-treats-Disease        : 1552 x 137           : 755         : 0.0036
hetio_CuG : Compound-upregulates-Gene      : 1552 x 20945         : 18756       : 5.8e-4
hetio_DaG : Disease-associates-Gene        : 137 x 20945          : 12623       : 0.0044
hetio_DdG : Disease-downregulates-Gene     : 137 x 20945          : 7623        : 0.0027
hetio_DlA : Disease-localizes-Anatomy      : 137 x 402            : 3602        : 0.0654
hetio_DpS : Disease-presents-Symptom       : 137 x 438            : 3357        : 0.0559
hetio_DrD : Disease-resembles-Disease      : 137 x 137            : 543         : 0.0289
hetio_DuG : Disease-upregulates-Gene       : 137 x 20945          : 7731        : 0.0027
hetio_GcG : Gene-covaries-Gene             : 20945 x 20945        : 61690       : 1.4e-4
hetio_GiG : Gene-interacs-Gene             : 20945 x 20945        : 147164      : 3.4e-4
hetio_GpBP: Gene-participates-Biological Process: 20945 x 11381   : 559504      : 0.0023
hetio_GpCC: Gene-participates-Cellular Component: 20945 x 1391    : 73566       : 0.0025
hetio_GpMF: Gene-participates-Molecular Function: 20945 x 2884    : 97222       : 0.0016
hetio_GpPW: Gene-participates-Pathway      : 20945 x 1822         : 84372       : 0.0022
hetio_GrG : Gene->regulates->gene (forward): 20945 x 20945        : 265672      : 6.1e-4
hetio_PCiC: Pharmacologic Class-includes-Compound: 345 x 1552     : 1029        : 0.0019
```
#### Similarity matrix(i,j)=float[0,1] for similarity between i and j
```
Name      : Edge type                      : Shape         : Method
hetio_CsC : Compound-similar-Compound      : 1552 x 1552   : JChem screenmd Tanimoto similarity of ECFP4
```
