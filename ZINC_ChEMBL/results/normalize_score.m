function N=normalize_score(normscore, P)
%normscore: array of normalized scores ordered by bin
%P: predicted raw score matrix
N=zeros(size(P,1),size(P,2));
for i=1:size(P,1)
 for j=1:size(P,2)
	  if P(i,j)<=1.5
		  raw=20*P(i,j);
		  f=fix(raw);
		  if (raw-f)==0
		   nscore=normscore(f);
		   N(i,j)=nscore;
		  else
		   nscore=normscore(f+1);
		   N(i,j)=nscore;
		  end
	  else
		   N(i,j)=P(i,j);	%raw score itself if greater than 1.5
	  end
 end
end
end
