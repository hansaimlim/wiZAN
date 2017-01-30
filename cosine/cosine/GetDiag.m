function [DA newA] = GetDiag(A,J)
    [m,n] = size(A);
    A(1:m+1:m*m) = 0;
    % sort rows in descending order
    [SH ind] = sort(A,2,'descend');
    TOP_HORIZ = zeros(m,m);
    
    for i = 1:m
        for j = 1:J
            % next largest element in row i is at position col
            col = ind(i,j);
            TOP_HORIZ(i,col) = A(i,col);
        end
    end
    DA = diag((sum(TOP_HORIZ,1) + sum(TOP_HORIZ,2)')/2);
    newA = (TOP_HORIZ + TOP_HORIZ')/2;
end
