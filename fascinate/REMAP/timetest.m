function [cputime, gputime]=timetest()

cputime=zeros(1,10);
gputime=zeros(1,10);
for i=1:10
 t=REMAP_chem_dis_vector();
 cputime(1,i)=t;
end
for i=1:10
 t=REMAP_chem_dis_gpu(dis_dis_sim_vector_norm);
 gputime(1,i)=t;
end
cmean=mean(cputime)
cstd=std(cputime)
gmean=mean(gputime)
gstd=std(gputime)

end
