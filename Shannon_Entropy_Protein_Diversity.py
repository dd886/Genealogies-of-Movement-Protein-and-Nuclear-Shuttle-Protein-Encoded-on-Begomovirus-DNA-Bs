
# coding: utf-8

# In[3]:


from Bio import SeqIO, AlignIO


# In[4]:


import math


# In[30]:


import pandas as pd
import trace


# In[6]:


import glob
import numpy as np
AAs = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']


# In[22]:


def shannon_info(site, AAs=AAs):
    
    temp_info = 0
    
    for aa in AAs:
        
        if aa in site:
            
            temp_info += site.count(aa)/len(site) * math.log2(site.count(aa)/len(site))
        
    return abs(temp_info)


# In[19]:


def log_20(x):
    try:
        return math.log(x,20)
    except:
        return 0

vlog_20 = np.vectorize(log_20)


# In[20]:


def vn_entropy(column,sub_matrix,matrix_label):

    column_diag = np.zeros((sub_matrix.shape[0],sub_matrix.shape[0]))
    
    for aa in matrix_label:
        column_diag[matrix_label.index(aa),matrix_label.index(aa)] = column.count(aa) / len(column)

    omega = column_diag * sub_matrix 
    
    entropy = - omega * vlog_20(omega)
    
    return entropy.trace()


# In[7]:


data_list = glob.glob("file name")


# In[8]:


LG_aa_order = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
B50_aa_order = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','J','Z','X','*']


# In[10]:


LG = pd.read_csv("LG.csv",header=0, index_col=0)


# In[11]:


B50 = pd.read_table("BLOSUM50.tab",header=None,index_col=None,names=B50_aa_order)
B50.set_index(pd.Series(B50_aa_order),inplace=True)


# In[12]:


max_freq = np.array([[0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.05,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                     [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]])


# In[31]:


#VNE value at Shannon Entropy Max: all amino acids at frequency 0.05
omega = B50 * max_freq
ent = -(omega * vlog_20(omega))


# In[32]:


ent.trace()
#move past this if attribute error returned


# In[35]:


print(ent)


# In[36]:


#Pairwise similarity maximum for two allele locus (C,E)
(1/LG.loc['C','E']) * 0.5


# In[40]:


#VNE maximum at two allele locus (C,W)
max_vne = np.array([[0.090909,0],[0,0.07878788]])
omega = max_vne * np.array([[0.5,0],[0,0.5]])
ent = -(omega * vlog_20(omega))


# In[41]:


#calculate the diversity measure
for currentfile in data_list:
    data = AlignIO.read(currentfile,'fasta')
    print("{} measures:".format(currentfile))    


    num_samples = len(data)
    ungapped_sites = 0

    shannon = 0
        
    #get diversity meaasure(s) for ungapped sites
    for i in range(len(data[0,:])):
    
        column = data[:,i]
        if '-' not in column:
        
            shannon += shannon_info(column)
    
    print("shannon: {:.4}".format(shannon/(i+1)))

