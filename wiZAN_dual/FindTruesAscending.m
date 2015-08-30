function [row_col_rank_score] = FindTruesAscending(A, TM) %A: predicted score matrix, TM: True Matrix (either true positive or true negative)
%work same as FindTrues.m, but find the rank in ascending order
%i.e. smaller values rank closer to the top (1st)
[i j k]=find(TM);
row_col_rank_score = zeros(length(i),4);
for n=1:length(i)
  row=i(n);
  col=j(n);
  predscore = A(row, col);
  drank = sum(A(row,:) <= A(row, col));
  row_col_rank_score(n,:) = [row, col, drank, predscore];
end 
end
