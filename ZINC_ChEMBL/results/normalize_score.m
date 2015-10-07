function normalize_score(normscore, P)
%normscore: array of normalized scores ordered by bin
%P: predicted raw score matrix
filehandles=[];
for i=1:size(normscore,2)
 nscore=normscore(i);
 filename=['./pairs_in_normScore_' num2str(nscore) '.txt'];
 fh=fopen(filename, 'a');
 filehandles=[filehandles, fh];
end

for i=1:size(P,1)
 for j=1:size(P,2)
  raw=20*P(i,j);
  f=fix(raw);
  if (raw-f)==0
   nscore=normscore(f);
   fh=filehandles(f);
   fpritf(fh, '%d, %d\n', [i, j]);
  else
   nscore=normscore(f+1);
   fh=filehandles(f+1);
   fpritf(fh, '%d, %d\n', [i, j]);
  end
 end
end
end
