function [RCRS] = FindTrues(A, TM) %A: predicted score matrix, TM: Test Matrix (either true positive or true negative)
%returns RCRS=[RowIndex, ColumnIndex, RowRank, PredScore]
[i j k]=find(TM);
RCRS = zeros(length(i),4);
for n=1:length(i)
  row=i(n);
  col=j(n);
  predscore = A(row, col);
  drank = sum(A(row,:) >= A(row, col));	%dense rank on the row
  RCRS(n,:) = [row, col, drank, predscore];
end 
end
