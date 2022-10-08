import pandas as pd
import numpy as np
# MC
df_MC=pd.read_table('EnhancerList_MC.txt', sep='\t')
df_MC = df_MC[df_MC['class'].str.contains('intergenic')]
df_MC1=df_MC[['chr','start', 'end', 'cellType', 'activity_base']]
df_MC2=df_MC[['chr','start', 'end', 'activity_base']]
df_MC1.to_csv('MC_enhancer_activity.csv')
df_MC2.to_csv('MC_enhancer_activity.bed', index=False,header=False,sep='\t')
# FL
df_FL=pd.read_table('EnhancerList_FL.txt', sep='\t')
df_FL = df_FL[df_FL['class'].str.contains('intergenic')]
df_FL1=df_FL[['chr','start', 'end', 'cellType', 'activity_base']]
df_FL2=df_FL[['chr','start', 'end', 'activity_base']]
df_FL1.to_csv('FL_enhancer_activity.csv')
df_FL2.to_csv('FL_enhancer_activity.bed', index=False,header=False,sep='\t')