#!/usr/bin/python

import os
import sys
import json
import numpy as np
import pandas as pd
from scipy.io import savemat
from scipy import sparse

def write_out(df,outfile,additional_id=None):
  idx=1
  with open(outfile,'w') as out:
    for dfidx, row in df.iterrows():
      if additional_id:
        dfid=str(row[0])
        addid=str(row[1])
        name=str(row[2])
        line=','.join((str(idx),dfid,addid,name))+'\n'
      else:
        dfid=str(row[0])
        name=str(row[1])
        line=','.join((str(idx),dfid,name))+'\n'
      out.write(line)
      idx+=1

def write_out_compound(df,outfile):
  idx=1
  with open(outfile, 'w') as out:
    for dfidx, row in df.iterrows():
      dfid=str(row[0])
      ikey=str(row[1])
      name=str(row[2])
      line=','.join((str(idx),dfid,ikey,name))+'\n'
      out.write(line)
      idx+=1

def write_out_gene(df,outfile):
  idx=1
  with open(outfile, 'w') as out:
    for dfidx, row in df.iterrows():
      dfid=str(row[0])
      name=str(row[1])
      desc=str(row[2])
      chro=str(row[3])
      line=','.join((str(idx),dfid,name,desc,chro))+'\n'
      out.write(line)
      idx+=1


def get_dict(df):
  #get dictionary of idx2id[index]=id and id2idx[id]=index
  idx=1
  idx2id={}
  id2idx={}
  for dfidx, row in df.iterrows():
    dfid=str(row[0])
    idx2id[idx]=dfid
    id2idx[dfid]=idx
    idx+=1
  return (idx2id, id2idx, len(id2idx))

def im1(indice):
  #index minus 1
  #reduce each index by 1 for matlab-python conversion
  idx=np.array(indice)-1
  return idx

def get_node(jsonfile,writeout=False):
  with open(jsonfile, "r") as inf:
    k=json.load(inf)
    df=pd.io.json.json_normalize(k)
    #extract nodes
    nodes=df['nodes']
    nodes_norm=pd.io.json.json_normalize(nodes[0])
    #anatomy list
    anatomy=nodes_norm.loc[nodes_norm['kind']=='Anatomy']
   # anatomy=anatomy[['data.mesh_id','identifier','name']]
    anatomy=anatomy[['identifier','name']]
    #bio process list
    BP=nodes_norm.loc[nodes_norm['kind']=='Biological Process']
    BP=BP[['identifier','name']]
    #cell component list
    CC=nodes_norm.loc[nodes_norm['kind']=='Cellular Component']
    CC=CC[['identifier','name']]
    #Compound list
    compound=nodes_norm.loc[nodes_norm['kind']=='Compound']
    compound=compound[['identifier','data.inchikey','name']]
    #disease list
    disease=nodes_norm.loc[nodes_norm['kind']=='Disease']
    disease=disease[['identifier','name']]
    #gene list
    gene=nodes_norm.loc[nodes_norm['kind']=='Gene']
    gene=gene[['identifier','name','data.description','data.chromosome']]
    #molec function list
    MF=nodes_norm.loc[nodes_norm['kind']=='Molecular Function']
    MF=MF[['identifier','name']]
    #pathway list
    pathway=nodes_norm.loc[nodes_norm['kind']=='Pathway']
    pathway=pathway[['identifier','name']]
    #pharmacologic class list
    PC=nodes_norm.loc[nodes_norm['kind']=='Pharmacologic Class']
   # PC=PC[['data.class_type','identifier','name']]
    PC=PC[['identifier','name']]
    #side effect list
    SE=nodes_norm.loc[nodes_norm['kind']=='Side Effect']
    SE=SE[['identifier','name']]
    #symptom list
    symp=nodes_norm.loc[nodes_norm['kind']=='Symptom']
    symp=symp[['identifier','name']]
    if writeout:
      anatfile='./list/hetio_anatomy_list.csv'
      write_out(anatomy,anatfile)
      bpfile='./list/hetio_bioprocess_list.csv'
      write_out(BP,bpfile)
      ccfile='./list/hetio_cellcomponent_list.csv'
      write_out(CC,ccfile)
      compfile='./list/hetio_compound_list.csv'
      write_out_compound(compound,compfile)
      disfile='./list/hetio_disease_list.csv'
      write_out(disease,disfile)
      genefile='./list/hetio_gene_list.csv'
      write_out_gene(gene,genefile)
      mffile='./list/hetio_molecfunction_list.csv'
      write_out(MF,mffile)
      pathfile='./list/hetio_pathway_list.csv'
      write_out(pathway,pathfile)
      pcfile='./list/hetio_pharmacoclass_list.csv'
      write_out(PC,pcfile,additional_id=True)
      sefile='./list/hetio_sideeffect_list.csv'
      write_out(SE,sefile)
      sympfile='./list/hetio_symptom_list.csv'
      write_out(symp,sympfile)
    #get dictionarys for ids and indices
    anat_idx2id, anat_id2idx, anat_max_idx = get_dict(anatomy)
    bp_idx2id, bp_id2idx, bp_max_idx = get_dict(BP)
    cc_idx2id, cc_id2idx, cc_max_idx = get_dict(CC)
    comp_idx2id, comp_id2idx, comp_max_idx = get_dict(compound)
    dis_idx2id, dis_id2idx, dis_max_idx = get_dict(disease)
    gene_idx2id, gene_id2idx, gene_max_idx = get_dict(gene)
    mf_idx2id, mf_id2idx, mf_max_idx = get_dict(MF)
    path_idx2id, path_id2idx, path_max_idx = get_dict(pathway)
    pc_idx2id, pc_id2idx, pc_max_idx = get_dict(PC)
    se_idx2id, se_id2idx, se_max_idx = get_dict(SE)
    symp_idx2id, symp_id2idx, symp_max_idx = get_dict(symp)
    #dictionary to return (contains idx2id and id2idx dicts
    DICTS={}
    DICTS['Anatomy']=(anat_idx2id, anat_id2idx, anat_max_idx)
    DICTS['Biological Process']=(bp_idx2id, bp_id2idx, bp_max_idx)
    DICTS['Cellular Component']=(cc_idx2id, cc_id2idx, cc_max_idx)
    DICTS['Compound']=(comp_idx2id, comp_id2idx, comp_max_idx)
    DICTS['Disease']=(dis_idx2id, dis_id2idx, dis_max_idx)
    DICTS['Gene']=(gene_idx2id, gene_id2idx, gene_max_idx)
    DICTS['Molecular Function']=(mf_idx2id, mf_id2idx, mf_max_idx)
    DICTS['Pathway']=(path_idx2id, path_id2idx, path_max_idx)
    DICTS['Pharmacologic Class']=(pc_idx2id, pc_id2idx, pc_max_idx)
    DICTS['Side Effect']=(se_idx2id, se_id2idx, se_max_idx)
    DICTS['Symptom']=(symp_idx2id, symp_id2idx, symp_max_idx)
    return DICTS
      
def get_edge(jsonfile, DICTS, writeout=False):
  #extract edges
  (anat_idx2id, anat_id2idx, anat_max_idx)=DICTS['Anatomy']
  (bp_idx2id, bp_id2idx, bp_max_idx)=DICTS['Biological Process']
  (cc_idx2id, cc_id2idx, cc_max_idx)=DICTS['Cellular Component']
  (comp_idx2id, comp_id2idx, comp_max_idx)=DICTS['Compound']
  (dis_idx2id, dis_id2idx, dis_max_idx)=DICTS['Disease']
  (gene_idx2id, gene_id2idx, gene_max_idx)=DICTS['Gene']
  (mf_idx2id, mf_id2idx, mf_max_idx)=DICTS['Molecular Function']
  (pw_idx2id, pw_id2idx, pw_max_idx)=DICTS['Pathway']
  (pc_idx2id, pc_id2idx, pc_max_idx)=DICTS['Pharmacologic Class']
  (se_idx2id, se_id2idx, se_max_idx)=DICTS['Side Effect']
  (symp_idx2id, symp_id2idx, symp_max_idx)=DICTS['Symptom']
  with open(jsonfile, "r") as inf:
    k=json.load(inf)
    df=pd.io.json.json_normalize(k)
    edges=df['edges']
    edges_norm=pd.io.json.json_normalize(edges[0])
    m,n=edges_norm.shape
    #source=Anatomy
    AdG_row=[];AdG_col=[]
    AeG_row=[];AeG_col=[]
    AuG_row=[];AuG_col=[]
    #source=Compound
    CbG_row=[];CbG_col=[]
    CdG_row=[];CdG_col=[]
    CuG_row=[];CuG_col=[]
    CcSE_row=[];CcSE_col=[]
    CpD_row=[];CpD_col=[]
    CtD_row=[];CtD_col=[]
    #source=Disease
    DaG_row=[];DaG_col=[]
    DdG_row=[];DdG_col=[]
    DuG_row=[];DuG_col=[]
    DlA_row=[];DlA_col=[]
    DpS_row=[];DpS_col=[]
    DrD_row=[];DrD_col=[]
    #source=Gene
    GcG_row=[];GcG_col=[]
    GiG_row=[];GiG_col=[]
    GpBP_row=[];GpBP_col=[]
    GpCC_row=[];GpCC_col=[]
    GpMF_row=[];GpMF_col=[]
    GpPW_row=[];GpPW_col=[]
    GrG_row=[];GrG_col=[]
    PCiC_row=[];PCiC_col=[]
    for i in range(m):
      source_type, source_id=edges_norm['source_id'][i]
      target_type, target_id=edges_norm['target_id'][i]
      direction=edges_norm['direction'][i]
      source_id=str(source_id);target_id=str(target_id)
      kind=edges_norm['kind'][i] #e.g. upregulate, downregulate, etc
      direction=edges_norm['direction'][i]
      if source_type=='Anatomy':
        #Anatomy-gene
        source_idx=anat_id2idx[source_id]
        target_idx=gene_id2idx[target_id]
        if kind=='downregulates':
          #AdG
          AdG_row.append(source_idx);AdG_col.append(target_idx)
        elif kind=='expresses':
          #AeG
          AeG_row.append(source_idx);AeG_col.append(target_idx)
        elif kind=='upregulates':
          #AuG
          AuG_row.append(source_idx);AuG_col.append(target_idx)
        else:
          #somethine else
          continue
      elif source_type=='Compound':
        #compound to gene/side effect/disease
        source_idx=comp_id2idx[source_id]
        if target_type=='Gene':
          #compound-gene
          target_idx=gene_id2idx[target_id]
          if kind=='binds':
            #CbG
            CbG_row.append(source_idx);CbG_col.append(target_idx)
          elif kind=='downregulates':
            #CdG
            CdG_row.append(source_idx);CdG_col.append(target_idx)
          elif kind=='upregulates':
            #CuG
            CuG_row.append(source_idx);CuG_col.append(target_idx)
          else:
            continue
        elif target_type=='Side Effect':
          #CcSE
          target_idx=se_id2idx[target_id]
          CcSE_row.append(source_idx);CcSE_col.append(target_idx)
        elif target_type=='Disease':
          target_idx=dis_id2idx[target_id]
          if kind=='palliates':
            #CpD
            CpD_row.append(source_idx);CpD_col.append(target_idx)
          elif kind=='treats':
            #CtD
            CtD_row.append(source_idx);CtD_col.append(target_idx)
          else:
            continue
      elif source_type=='Disease':
        #disease to gene/anatomy/symptom/disease
        source_idx=dis_id2idx[source_id]
        if kind=='associates':
          #DaG
          target_idx=gene_id2idx[target_id]
          DaG_row.append(source_idx);DaG_col.append(target_idx)
        elif kind=='upregulates':
          #DuG
          target_idx=gene_id2idx[target_id]
          DuG_row.append(source_idx);DuG_col.append(target_idx)
        elif kind=='downregulates':
          #DdG
          target_idx=gene_id2idx[target_id]
          DdG_row.append(source_idx);DdG_col.append(target_idx)
        elif kind=='localizes':
          #DlA
          target_idx=anat_id2idx[target_id]
          DlA_row.append(source_idx);DlA_col.append(target_idx)
        elif kind=='presents':
          #DpS
          target_idx=symp_id2idx[target_id]
          DpS_row.append(source_idx);DpS_col.append(target_idx)
        elif kind=='resembles':
          #DrD
          target_idx=dis_id2idx[target_id]
          DrD_row.append(source_idx);DrD_col.append(target_idx)
        else:
          continue
      elif source_type=='Gene':
        #gene to gene/biological process/molecular function/pathway
        source_idx=gene_id2idx[source_id]
        if target_type=='Gene':
          target_idx=gene_id2idx[target_id]
          if kind=='covaries':
            #GcG
            GcG_row.append(source_idx);GcG_col.append(target_idx)
          elif kind=='interacts':
            #GiG
            GiG_row.append(source_idx);GiG_col.append(target_idx)
          elif kind=='regulates':
            #Gr>G ###directional
            GrG_row.append(source_idx);GrG_col.append(target_idx)
          else:
            continue
        elif target_type=='Biological Process':
          #GpBP
          target_idx=bp_id2idx[target_id]
          GpBP_row.append(source_idx);GpBP_col.append(target_idx)
        elif target_type=='Cellular Component':
          #GpCC
          target_idx=cc_id2idx[target_id]
          GpCC_row.append(source_idx);GpCC_col.append(target_idx)
        elif target_type=='Molecular Function':
          #GpMF
          target_idx=mf_id2idx[target_id]
          GpMF_row.append(source_idx);GpMF_col.append(target_idx)
        elif target_type=='Pathway':
          #GpPW
          target_idx=pw_id2idx[target_id]
          GpPW_row.append(source_idx);GpPW_col.append(target_idx)
        else:
          continue
      elif source_type=='Pharmacologic Class':
        source_idx=pc_id2idx[source_id]
        if target_type=='Compound':
          #PCiC
          target_idx=comp_id2idx[target_id]
          PCiC_row.append(source_idx);PCiC_col.append(target_idx)
        else:
          continue
      else:
        continue
    #matlab-type indice to python-type indice
    AdG_row=im1(AdG_row);AdG_col=im1(AdG_col)
    AeG_row=im1(AeG_row);AeG_col=im1(AeG_col)
    AuG_row=im1(AuG_row);AuG_col=im1(AuG_col)
    #source=Compound
    CbG_row=im1(CbG_row);CbG_col=im1(CbG_col)
    CdG_row=im1(CdG_row);CdG_col=im1(CdG_col)
    CuG_row=im1(CuG_row);CuG_col=im1(CuG_col)
    CcSE_row=im1(CcSE_row);CcSE_col=im1(CcSE_col)
    CpD_row=im1(CpD_row);CpD_col=im1(CpD_col)
    CtD_row=im1(CtD_row);CtD_col=im1(CtD_col)
    #source=Disease
    DaG_row=im1(DaG_row);DaG_col=im1(DaG_col)
    DdG_row=im1(DdG_row);DdG_col=im1(DdG_col)
    DuG_row=im1(DuG_row);DuG_col=im1(DuG_col)
    DlA_row=im1(DlA_row);DlA_col=im1(DlA_col)
    DpS_row=im1(DpS_row);DpS_col=im1(DpS_col)
    DrD_row=im1(DrD_row);DrD_col=im1(DrD_col)
    #source=Gene
    GcG_row=im1(GcG_row);GcG_col=im1(GcG_col)
    GiG_row=im1(GiG_row);GiG_col=im1(GiG_col)
    GpBP_row=im1(GpBP_row);GpBP_col=im1(GpBP_col)
    GpCC_row=im1(GpCC_row);GpCC_col=im1(GpCC_col)
    GpMF_row=im1(GpMF_row);GpMF_col=im1(GpMF_col)
    GpPW_row=im1(GpPW_row);GpPW_col=im1(GpPW_col)
    GrG_row=im1(GrG_row);GrG_col=im1(GrG_col)
    PCiC_row=im1(PCiC_row);PCiC_col=im1(PCiC_col)
    
    AdG=sparse.coo_matrix(([1]*len(AdG_row),(AdG_row,AdG_col)), shape=(anat_max_idx,gene_max_idx),dtype=np.int16)
    AeG=sparse.coo_matrix(([1]*len(AeG_row),(AeG_row,AeG_col)), shape=(anat_max_idx,gene_max_idx),dtype=np.int16)
    AuG=sparse.coo_matrix(([1]*len(AuG_row),(AuG_row,AuG_col)), shape=(anat_max_idx,gene_max_idx),dtype=np.int16)
    CbG=sparse.coo_matrix(([1]*len(CbG_row),(CbG_row,CbG_col)), shape=(comp_max_idx,gene_max_idx),dtype=np.int16)
    CdG=sparse.coo_matrix(([1]*len(CdG_row),(CdG_row,CdG_col)), shape=(comp_max_idx,gene_max_idx),dtype=np.int16)
    CuG=sparse.coo_matrix(([1]*len(CuG_row),(CuG_row,CuG_col)), shape=(comp_max_idx,gene_max_idx),dtype=np.int16)
    CcSE=sparse.coo_matrix(([1]*len(CcSE_row),(CcSE_row,CcSE_col)), shape=(comp_max_idx,se_max_idx),dtype=np.int16)
    CpD=sparse.coo_matrix(([1]*len(CpD_row),(CpD_row,CpD_col)), shape=(comp_max_idx,dis_max_idx),dtype=np.int16)
    CtD=sparse.coo_matrix(([1]*len(CtD_row),(CtD_row,CtD_col)), shape=(comp_max_idx,dis_max_idx),dtype=np.int16)
    DaG=sparse.coo_matrix(([1]*len(DaG_row),(DaG_row,DaG_col)), shape=(dis_max_idx,gene_max_idx),dtype=np.int16)
    DdG=sparse.coo_matrix(([1]*len(DdG_row),(DdG_row,DdG_col)), shape=(dis_max_idx,gene_max_idx),dtype=np.int16)
    DuG=sparse.coo_matrix(([1]*len(DuG_row),(DuG_row,DuG_col)), shape=(dis_max_idx,gene_max_idx),dtype=np.int16)
    DlA=sparse.coo_matrix(([1]*len(DlA_row),(DlA_row,DlA_col)), shape=(dis_max_idx,anat_max_idx),dtype=np.int16)
    DpS=sparse.coo_matrix(([1]*len(DpS_row),(DpS_row,DpS_col)), shape=(dis_max_idx,symp_max_idx),dtype=np.int16)
    DrD=sparse.coo_matrix(([1]*len(DrD_row),(DrD_row,DrD_col)), shape=(dis_max_idx,dis_max_idx),dtype=np.int16)
    GcG=sparse.coo_matrix(([1]*len(GcG_row),(GcG_row,GcG_col)), shape=(gene_max_idx,gene_max_idx),dtype=np.int16)
    GiG=sparse.coo_matrix(([1]*len(GiG_row),(GiG_row,GiG_col)), shape=(gene_max_idx,gene_max_idx),dtype=np.int16)
    GpBP=sparse.coo_matrix(([1]*len(GpBP_row),(GpBP_row,GpBP_col)), shape=(gene_max_idx,bp_max_idx),dtype=np.int16)
    GpCC=sparse.coo_matrix(([1]*len(GpCC_row),(GpCC_row,GpCC_col)), shape=(gene_max_idx,cc_max_idx),dtype=np.int16)
    GpMF=sparse.coo_matrix(([1]*len(GpMF_row),(GpMF_row,GpMF_col)), shape=(gene_max_idx,mf_max_idx),dtype=np.int16)
    GpPW=sparse.coo_matrix(([1]*len(GpPW_row),(GpPW_row,GpPW_col)), shape=(gene_max_idx,pw_max_idx),dtype=np.int16)
    GrG=sparse.coo_matrix(([1]*len(GrG_row),(GrG_row,GrG_col)), shape=(gene_max_idx,gene_max_idx),dtype=np.int16)
    PCiC=sparse.coo_matrix(([1]*len(PCiC_row),(PCiC_row,PCiC_col)), shape=(pc_max_idx,comp_max_idx),dtype=np.int16)
    if writeout:
      savemat('./mat/hetio_AdG', {'hetio_AdG':AdG})
      savemat('./mat/hetio_AeG', {'hetio_AeG':AeG})
      savemat('./mat/hetio_AuG', {'hetio_AuG':AuG})
      savemat('./mat/hetio_CbG', {'hetio_CbG':CbG})
      savemat('./mat/hetio_CdG', {'hetio_CdG':CdG})
      savemat('./mat/hetio_CuG', {'hetio_CuG':CuG})
      savemat('./mat/hetio_CcSE', {'hetio_CcSE':CcSE})
      savemat('./mat/hetio_CpD', {'hetio_CpD':CpD})
      savemat('./mat/hetio_CtD', {'hetio_CtD':CtD})

      savemat('./mat/hetio_DaG', {'hetio_DaG':DaG})
      savemat('./mat/hetio_DdG', {'hetio_DdG':DdG})
      savemat('./mat/hetio_DuG', {'hetio_DuG':DuG})
      savemat('./mat/hetio_DlA', {'hetio_DlA':DlA})
      savemat('./mat/hetio_DpS', {'hetio_DpS':DpS})
      savemat('./mat/hetio_DrD', {'hetio_DrD':DrD})

      savemat('./mat/hetio_GcG', {'hetio_GcG':GcG})
      savemat('./mat/hetio_GiG', {'hetio_GiG':GiG})
      savemat('./mat/hetio_GpBP', {'hetio_GpBP':GpBP})
      savemat('./mat/hetio_GpCC', {'hetio_GpCC':GpCC})
      savemat('./mat/hetio_GpMF', {'hetio_GpMF':GpMF})
      savemat('./mat/hetio_GpPW', {'hetio_GpPW':GpPW})
      savemat('./mat/hetio_GrG', {'hetio_GrG':GrG})
      savemat('./mat/hetio_PCiC', {'hetio_PCiC':PCiC})

if __name__ == '__main__':
  jsonfile='./hetionet-v1.0.json'
  node_dicts=get_node(jsonfile)
  get_edge(jsonfile, node_dicts, writeout=True)
  
  
