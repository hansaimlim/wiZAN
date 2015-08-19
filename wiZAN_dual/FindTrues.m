function [row_col_rank_score] = FindTrues(A, TM) %A: predicted score matrix, TM: True Matrix (either true positive or true negative)
si = (size(A));
row_col_rank_score = zeros(sum(TM(:)),4);
rcrs_count = 1;
[i j k]=find(TM);
for n=1:length(i)
  row=i(n);
  col=j(n);
  predscore = A(row, col);
  drank = sum(A(row,:) >= A(row, col));
  row_col_rank_score(rcrs_count,:) = [row, col, drank, predscore];
  rcrs_count = rcrs_count + 1;
end 
end
