function [] = ZincTest()
    % add path to zinc files here
    % addpath('');
    load ('Zinc_R.dat');
    R = spconvert(Zinc_R);

    load ('Zinc_M.dat');
    M = spconvert(Zinc_M);
    
    load ('Zinc_N.dat');
    N = spconvert(Zinc_N);
    
    
    [m,n] = size(R);

    iter = 600; 
    J = 3;
    QUICK = 0;
    M_cut = 0.3;
    N_cut = 0.5;
    
    rnk=150; wp=1.0; lR=0.1; lM=1.0; lN=10.0; 

    for TEST=0:4
        testfile = strcat(num2str(TEST),'.test');
        delimiterIn = ' ';
        headerlinesIn = 0;
        T = importdata(testfile,delimiterIn,headerlinesIn);
    
        outfile = ['zinc_', num2str(TEST), '.out'];
        fid=fopen(outfile,'w');
    
        [TPR time] = TprBench(R,M,N,T,J,lR,lM,lN,iter,rnk,wp,M_cut,N_cut,QUICK);
        fprintf(fid,'TPR:%f rank:%d WP:%f lR:%f lM:%f lN:%f time:%f\n',TPR,rnk,wp,lR,lM,lN,time);
        fclose(fid);
    end

end
