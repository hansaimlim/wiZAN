This package enables you to run either cross validation on Yamanishi test sets or 
TPR analysis on ZINC test sets (see the manuscript for details). 


Extract:
%gunzip code.tar.gz
%tar -xvf code.tar
%cd code

To run the programs from the command line, enter, for instance:

% nohup matlab -nodesktop -nodisplay -nosplash -r "Run();exit" > stout >& sterr &

For Yamanishi, use Run() function; for ZINC uze ZincTest().

The path have to be set in these MATLAB scripts to point to correct folders containing 
Yamanishi or ZINC data.

## For actual prediction, please use COSINE_Predict()
