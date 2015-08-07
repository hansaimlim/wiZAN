function wiZAN_predscores(u_mat_csv, v_mat_csv, T, outfile)
%u_mat_csv is filename for the low rank matrix U in csv file
%U and V are low rank matrices; Pred=U*V'
%T is the matrix of interesting pairs (e.g. ambiguous pairs, or unknown pairs)
U=csvread(u_mat_csv);
V=csvread(v_mat_csv);
Pred=U*V';
clear U V;

Pred_t=Pred.*T;
clear Pred T;

[i, j, val]=find(Pred_t);
csvwrite(outfile, [val]);
clear
end
