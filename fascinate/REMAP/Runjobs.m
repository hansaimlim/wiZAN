function Runjobs()
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\chem_chem_sim.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\prot_prot_sim.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_blank.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_lin_norm.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_path_norm.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_resnik_norm.mat');
load('C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\data\similarity\dis_dis_sim_vector_norm.mat');
REMAP_chem_dis_blank();
REMAP_chem_dis(dis_dis_sim, 'blank', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_cd_blank.txt');
REMAP_chem_dis(dis_dis_sim_lin_norm, 'lin_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_cd_lin_norm.txt');
REMAP_chem_dis(dis_dis_sim_path_norm, 'path_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_cd_path_norm.txt');
REMAP_chem_dis(dis_dis_sim_resnik_norm, 'resnik_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_cd_resnik_norm.txt');
REMAP_chem_dis(dis_dis_sim_vector_norm, 'vector_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_cd_vector_norm.txt');

REMAP_gene_dis(dis_dis_sim, 'blank', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_gd_blank.txt');
REMAP_gene_dis(dis_dis_sim_lin_norm, 'lin_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_gd_lin_norm.txt');
REMAP_gene_dis(dis_dis_sim_path_norm, 'path_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_gd_path_norm.txt');
REMAP_gene_dis(dis_dis_sim_resnik_norm, 'resnik_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_gd_resnik_norm.txt');
REMAP_gene_dis(dis_dis_sim_vector_norm, 'vector_norm', 'C:\Users\Hansaim\Documents\GitHub\wiZAN\fascinate\REMAP\dis_sim\TPR_gd_vector_norm.txt');
end