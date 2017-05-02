function simMat=TaniSim(featMat)
%use to get Tanimoto similarity matrix from feature matrix
%for example, PaDeL-generated chemical feature csv file contains
%chemical features for a chemical on each row
%the number of rows is the matrix size (m by m)
%the output is symmetric similarity matrix with 1 on diagonals
%user may want to filter by cutoff score (e.g. sim > 0.5)
%no NaN is allowed


%features on columns
[m,n]=size(featMat);
simMat=zeros(m,m);

for i=1:m
    vec1=featMat(i,:);
    for j=1:m
        if j>i
            vec2=featMat(j,:);
            sumvec=vec1+vec2;
            simMat(i,j)=(vec1*vec2')/sum(sumvec>0);
        end 
    end
end
simMat(isnan(simMat))=0;
simMat=simMat+simMat'+eye(m);

end