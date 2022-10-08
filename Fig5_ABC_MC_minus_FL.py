import pandas as pd
df_MC_vs_FL=pd.read_table('MC_vs_FL_ABC_score.bed ', sep='\t')
new_col = ['MC_chr', 'MC_start','MC_end','ABC_MC','FL_chr', 'FL','FL_end','ABC_FL','overlap_length']
df_MC_vs_FL.columns=new_col
df_MC_vs_FL['ABC_MC-ABC_FL']=df_MC_vs_FL['ABC_MC']-(df_MC_vs_FL['ABC_FL'])
df_MC_vs_FL.to_csv('MC_vs_FL_ABC_score_differences.csv',index=False)