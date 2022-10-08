# Calculate MC_ABC_score
import pandas as pd
df_MC_ABC=pd.read_table('MC_enhancer_contact.bed', sep='\t')
new_col = ['Enh_chr', 'Enh_start','Enh_end','activity_MC','Contact_chr', 'Contact_start','Contact_end','contact_MC','overlap_length']
df_MC_ABC.columns=new_col
df_MC_ABC = df_MC_ABC.drop_duplicates(subset=['Enh_start'], keep='first')
df_MC_ABC['MC_AxC']=df_MC_ABC['contact_MC'].mul(df_MC_ABC['activity_MC'])
df_MC_ABC = df_MC_ABC.dropna()
df_MC_ABC['MC_ABC']=df_MC_ABC['MC_AxC'] / sum(df_MC_ABC['MC_AxC'])
df_MC_ABC.to_csv('MC_ABC_score.csv')
df_MC_ABC=df_MC_ABC[['Enh_chr', 'Enh_start','Enh_end','MC_ABC']]
df_MC_ABC.to_csv('MC_ABC_score.bed', header = False, index =  False, sep = '\t')
# Calculate FL_ABC_score
import pandas as pd
df_FL_ABC=pd.read_table('FL_enhancer_contact.bed', sep='\t')
new_col = ['Enh_chr', 'Enh_start','Enh_end','activity_FL','Contact_chr', 'Contact_start','Contact_end','contact_FL','overlap_length']
df_FL_ABC.columns=new_col
df_FL_ABC = df_FL_ABC.drop_duplicates(subset=['Enh_start'], keep='first')
df_FL_ABC['FL_AxC']=df_FL_ABC['contact_FL'].mul(df_FL_ABC['activity_FL'])
df_FL_ABC = df_FL_ABC.dropna()
df_FL_ABC['FL_ABC']=df_FL_ABC['FL_AxC'] / sum(df_FL_ABC['FL_AxC'])
df_FL_ABC.to_csv('FL_ABC_score.csv')
df_FL_ABC=df_FL_ABC[['Enh_chr', 'Enh_start','Enh_end','FL_ABC']]
df_FL_ABC.to_csv('FL_ABC_score.bed', header = False, index =  False, sep = '\t')