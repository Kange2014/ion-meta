#!/usr/bin/env python

# -*- coding:utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import click

def arg_png(out_png, arg_pd):
    #1 药物类型饼图
    drug_class=arg_pd['predicted_ARG-class']
    drug_class_type=drug_class.value_counts()
    drug_class_df=pd.DataFrame(drug_class_type)
    if arg_pd.empty:
        autopct = None
    else:
        autopct = '%3.1f%%'
    plt.pie(drug_class_df['predicted_ARG-class'],labels=drug_class_df.index, autopct = autopct)
    plt.savefig(out_png)
    #return out_png

#得到数据库中每个ARG的平均长度
def read_amr_len(ARG_length):
    d_arg_len={}
    d_arg_avg_len={}
    with open (ARG_length,'r') as f:
        for line in f:
            ARG, length = line.rstrip().split()
            ARG = ARG.split("|")[-1]
            ARG=ARG.upper()
            d_arg_len.setdefault(ARG,[]).append(int(length))
            
    for k,v in d_arg_len.items():
        avg = int(np.mean(v))
        d_arg_avg_len[k]=avg
    return d_arg_avg_len

#去除覆盖区间的overlap区域,计算ARG的cover region
def cover_region_depth(num_list):
    '''
    输入：线段起点、终点
    输出：总的覆盖长度
    '''
    start=10000
    end=0
    for one_list in num_list:
        if one_list[0]<start:
            start=one_list[0]
        if one_list[1]>end:
            end=one_list[1]
    #print start, end
    flag=['false']*end
    for i in range(len(num_list)):
        for j in range(num_list[i][0], num_list[i][1]):
            flag[j]=True
    #print flag
    return flag.count(True)

def arg_table(deeparg_file, arg_pd, d_arg_avg_len):
    #2 drug class ARG的reads数
    drug_class_arg=arg_pd.loc[:, ['predicted_ARG-class','#ARG','counts']]
    arg_count=drug_class_arg.groupby(['predicted_ARG-class','#ARG'],as_index=False).sum()

    #3 计算cover、深度、相似度
    d_arg_cover={}
    d_arg_identity={}
    with open (deeparg_file,'r') as f:
        next(f)  # skip the first line.
        for line in f:    
            info=line.rstrip().split()
            #temp=[]
            temp=[int(info[1]),int(info[2])]
            #d_arg_reads[info[0]]+=1
            d_arg_cover.setdefault(info[0],[]).append(temp)
            d_arg_identity.setdefault(info[0],[]).append(float(info[7]))
    #d_arg_avg_len=read_amr_len(test2)

    temp_list = []
    for arg in d_arg_cover.keys():
        cover_len=cover_region_depth(d_arg_cover[arg])
        cover_p=format(cover_len*100/d_arg_avg_len[arg],'.1f')
        id=max(d_arg_identity[arg])
        #cover_p=cover_len/read_amr_len(arg)
        #print("%s\t%d\t%.2f%%\t%.2f%%" % (arg, cover_region_depth(d_arg_cover[arg]),cover_p,id))
        temp_list.append([arg, cover_region_depth(d_arg_cover[arg]),d_arg_avg_len[arg],cover_p,id])
    cover_id_pd=pd.DataFrame(temp_list,columns=['#ARG', 'covered_length','gene_length','perc_covered', 'perc_identity'])
    final=pd.merge(arg_count, cover_id_pd, on='#ARG')
    final=final.rename(columns={'counts':'num_reads','#ARG':'ARG', 'gene_length': 'ARG_length'})
    return final

@click.command()
@click.option('--deeparg_file', help='the deeparg results as: XXX.mapping.ARG')
@click.option('--out_table', help='The output ARG table.')
@click.option('--out_png', help='The output ARG drug class png.')
@click.option('--features_gene_length', default="./data/database/v2/features.gene.length", help='The database ARG length file.')

def main(features_gene_length, deeparg_file, out_png, out_table):
    #features_gene_length="./data/database/v2/features.gene.length"
    d_arg_avg_len=read_amr_len(features_gene_length)
    arg_pd=pd.read_table(deeparg_file, header=0)
    arg_png(out_png, arg_pd)
    final=arg_table(deeparg_file, arg_pd, d_arg_avg_len)
    #xls=open(out_table,'w')
    #xls.write(final)
    final.to_csv(out_table,sep="\t",index=False)

if __name__=="__main__":
    main()