# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 11:27:55 2024

@author: Student
"""


# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 15:07:22 2024

@author: Student
"""
import random
import pandas as pd
import numpy as np
import os

'''*************** SETTINGS ***************'''

version='lineage exp s=e+0.05'

network_size=300        
no_connects=3           
time=300                
seed_species=6         

'''loop settings'''
runs=10

rs='s=e+0.05'

spec_count=0
ext_count=0



'''******************* FUNCTIONS *********************'''

# SPAWNING THE NETWORK
def init_network(network_size):
    nodes=list(range(network_size))
    targets=[]
    for n in nodes:
        possible_connections = list(range(network_size))
        possible_connections.remove(n)
        options=random.sample(possible_connections,no_connects)
        targets.append(options)
    for no,ta in zip(nodes,targets):
        for t in ta:
            if no not in targets[t]:
                targets[t].append(no)
   
    #sp=[random.choice([True,False]) for _ in range(network_size)]
    
    sp=[False for _ in range(network_size)]
    seeded = 0
    while seeded < seed_species:
        index=random.choice(range(network_size))
        while sp[index] == True:
            index=random.choice(range(network_size))
        sp[index]=True
        seeded = seeded + 1
    
    counter=0
    tt=[[] for p in sp]
    
    
    lsem=[(random.random() / 20) for _ in range(len(sp))]
    lssm=[]
    for i,e in enumerate(lsem):
        s=e+0.05
        lssm.append(s)

    for o,(t,v) in enumerate(zip(lssm,lsem)):
        if sp[o]==False:
            lssm[o]=0
            lsem[o]=0
        if sp[o]==True:
            counter=counter+1
            tt[o]=counter
    ln=['' for _ in range(network_size)]
    data={'node id':nodes,'target':targets,'species':sp,'lineage specific spec. rate mod':lssm,'lineage specific ext. rate mod':lsem ,'lineage':ln,'species id':tt}
    df=pd.DataFrame(data)
    dfex=df.explode('target')
    
   
            
    #print(df['lineage specific spec. rate mod'])
    return(df,dfex,counter)



    


def get_neighbours(df,node):            # get list of connections per node
    pairs=[(row['node id'], row['target']) for index, row in df.iterrows()]
    neighbourhoods={}
    for no,ta in pairs:
        if no not in neighbourhoods:
            neighbourhoods[no] = []
        neighbourhoods[no].append(ta)
    neighbourhoods = [[no] + ta for no, ta in neighbourhoods.items()]
    return(neighbourhoods[node])


def speciation(df,dcopy,node,rate,counter,spec_count): # Looks correct
    p=random.uniform(0,1)                   # random number
    to_speciate=random.choice(neighbours)  
    check=True
  
                
    if p<rate and (dcopy.loc[to_speciate,'species']).all()==False:
        if check==True:
            dcopy.loc[to_speciate,'species']=True
            spec_count=spec_count+1
            lineage_list=[df.loc[node,'species id']]# ENF - I've made a change here, so that lineage tracks species ID history rather than node history
            if isinstance(lineage_list, list):
                lineage = ','.join(map(str, lineage_list))
                
            lineage = f"{lineage}"
            dcopy.at[to_speciate, 'lineage'] = lineage
            counter=counter+1
            df.loc[to_speciate,'species id']=counter
            dcopy.loc[to_speciate,'lineage specific spec. rate mod']=df.loc[node,'lineage specific spec. rate mod']
            dcopy.loc[to_speciate,'lineage specific ext. rate mod']=df.loc[node,'lineage specific ext. rate mod']
        else: 
            return(dcopy,counter,spec_count)
                
    return(dcopy,counter,spec_count)


def extinction(dcopy,node,rate,ext_count): # Looks correct
    extinct=None                # True/False label for each node
    check=True

    p=random.uniform(0,1)
    if p<rate:
        dcopy.loc[node,'species']=False
        ext_count=ext_count+1
        dcopy.loc[node,'lineage']=''
        dcopy.loc[node,'lineage specific spec. rate mod']=0
        dcopy.loc[node,'lineage specific ext. rate mod']=0
        df.loc[node,'species id']=''
       
        extinct=True
    if p>rate:
        extinct=False
    
   
    if check==False:                # if any requirement not met, go extinct
        dcopy.loc[node,'species']=False
        ext_count=ext_count+1
        dcopy.loc[node,'lineage']=''
        dcopy.loc[node,'lineage specific spec. rate mod']=0
        dcopy.loc[node,'lineage specific ext. rate mod']=0
        df.loc[node,'species tracker']=''
        extinct=True
    
    return(dcopy,extinct,ext_count)
 

'''***************** CALLING FUNCTIONS ****************'''

timeline=[]             

''' LOOP '''

folder_name='{}'.format(version)

os.makedirs(folder_name, exist_ok=True)

for r in range(8,runs):
    empty=0
    
    df,dfex,counter=init_network(network_size)
    
    file_name='{}.{}.csv'.format(version,r)
    file_path=os.path.join(folder_name, file_name)
    
    settings={'RS ':rs}
    settings_tab=pd.DataFrame(settings.items(), columns=['',''])
    
    settings_tab.to_csv(file_path,mode='w', index=False)

    with open(file_path, 'a') as f:
        f.write('\n# Original network\n')
    
    df.to_csv(file_path,mode='a',index=False)
    timeline.append(dfex)
    
    species_network = df[df['species'] == True]
    species_count = species_network.shape[0]

    print('TIME = ',0)
    time_title=('TIME = ',0)
    time_stats={'@ TIME':0,'Number of speciation events':spec_count,'Number of extinction events':ext_count,'Number of species':species_count}
    time_tab=pd.DataFrame(time_stats.items(),columns=['',''])
    time_tab.to_csv(file_path,mode='a', index=False)
    with open(file_path, 'a') as f:
        f.write('# Species Network Data\n')

    species_network.to_csv(file_path,mode='a', index=False)


    for t in range(time):
        spec_count=0
        ext_count=0
        dfex=df.explode('target')
        dcopy=dfex.copy()   
            
        for node,row in df.iterrows():
            
            lineage_sm=df.iloc[node,3]                
            lineage_em=df.iloc[node,4]
           
            species_status=df.loc[node,'species']
            
            neighbours=get_neighbours(dfex,node)
            
            if species_status==True:
                dcopy,counter,spec_count=speciation(df,dcopy,node,lineage_sm,counter,spec_count)
                dcopy,extinct,ext_count=extinction(dcopy,node,lineage_em,ext_count)
             
        # update network with changes during this time step:
        small_df=dcopy.groupby('node id')['target'].apply(list).reset_index()
        d2=dcopy.groupby('node id')['species'].any().reset_index()
        df.set_index('node id', inplace=True)
        d2.set_index('node id', inplace=True)
        df['species'] = d2['species']
        d2=dcopy.groupby('node id')['lineage'].first().reset_index()
        df['lineage']=d2['lineage']
        d2=dcopy.groupby('node id')['lineage specific spec. rate mod'].first().reset_index()
        df['lineage specific spec. rate mod']=d2['lineage specific spec. rate mod']
        d2=dcopy.groupby('node id')['lineage specific ext. rate mod'].first().reset_index()
        df['lineage specific ext. rate mod']=d2['lineage specific ext. rate mod']
        df.reset_index(inplace=True)
        species_network = df[df['species'] == True]
     
        species_count = species_network.shape[0]
        
        t=t+1

        print('TIME = ',t)
        time_title=('TIME = ',t)
        
        time_stats={'@ TIME':t,'Number of speciation events':spec_count,'Number of extinction events':ext_count,'Number of species':species_count}
        time_tab=pd.DataFrame(time_stats.items(),columns=['',''])
        
        time_tab.to_csv(file_path,mode='a', index=False)
       
        with open(file_path, 'a') as f:
            f.write('# Species Network Data\n')
        
        species_network.to_csv(file_path,mode='a', index=False)
        timeline.append(dcopy) 
        if species_count==0:
            empty=empty+1
        if empty==2:
            break
