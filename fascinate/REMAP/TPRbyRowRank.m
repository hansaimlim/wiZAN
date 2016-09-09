function [TPRbyRow] = TPRbyRowRank(RCRS, maxRank)
%to get TPR by cutoff Row rank
%output contains vector of dimension (maxRank, 2)
%Each row of RCRS contains one true positive association: [RowIndex,ColumnIndex,RowRank,PredScore]
%RCRS input must contain row rank on the 3rd column
%RCRS input must contain ONLY test pairs, since the number of lines will be used as condition positive
[cp, col] = size(RCRS); %cp: condition positive, col: number of columns
TPRbyRow = zeros(maxRank, 2);
row_rank = RCRS(:, 3);    %Rank by Rows on 3rd column
rank = 1;       %cutoff rank
while rank <= maxRank
        tpr = sum(row_rank<=rank)/cp;
        TPRbyRow(rank, :) = [rank, tpr];        %update result row
        rank = rank + 1;
end
end
