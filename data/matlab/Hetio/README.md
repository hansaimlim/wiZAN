### Data matrices for Hetionet: https://github.com/dhimmel/hetionet.

### The matlab matrices can be imported from python.
```py
from scipy.io import savemat, loadmat
x=loadmat('hetio_AdG.mat')
```
### Matrices information
--- matrix(i,j)=1 if ith row and jth column is positive
```
Name      : Edge type                      : Shape
hetio_AdG : Anatomy-downregulates-Gene     : 402 x 20945
hetio_AeG : Anatomy-expresses-Gene         : 402 x 20945
hetio_AuG : Anatomy-upregulates-Gene       : 402 x 20945
hetio_CbG : Compound-binds-Gene            : 1552 x 20945
hetio_CcSE: Compound-causes-Side Effect    : 1552 x 5734
hetio_CdG : Compound-downregulates-Gene    : 1552 x 20945
hetio_CpD : Compound-palliates-Disease     : 1552 x 137
hetio_CtD : Compound-treats-Disease        : 1552 x 137
hetio_CuG : Compound-upregulates-Gene      : 1552 x 20945
hetio_DaG : Disease-associates-Gene        : 137 x 20945
hetio_DdG : Disease-downregulates-Gene     : 137 x 20945
hetio_DlA : Disease-localizes-Anatomy      : 137 x 402
hetio_DpS : Disease-presents-Symptom       : 137 x 438
hetio_DrD : Disease-resembles-Disease      : 137 x 137
hetio_DuG : Disease-upregulates-Gene       : 137 x 20945
hetio_GcG : Gene-covaries-Gene             : 20945 x 20945
hetio_GiG : Gene-interacs-Gene             : 20945 x 20945
hetio_GpBP: Gene-participates-Biological Process: 20945 x 11381
hetio_GpCC: Gene-participates-Cellular Component: 20945 x 1391
hetio_GpMF: Gene-participates-Molecular Function: 20945 x 2884
hetio_GpPW: Gene-participates-Pathway      : 20945 x 1822
hetio_GrG : Gene->regulates->gene (forward): 20945 x 20945
hetio_PCiC: Pharmacologic Class-includes-Compound: 345 x 1552
```
