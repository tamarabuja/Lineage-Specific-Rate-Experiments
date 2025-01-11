# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 14:22:42 2024

@author: Student
"""




import random
import pandas as pd
#import numpy as np
import networkx as nx
#from pyvis.network import Network
#from IPython.display import HTML
import pathpy as pp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
#from functools import reduce
#import copy
import string
import os


#TEST VER
#random.seed(6)
# remove this to introduce randomness


'''*************** SETTINGS ***************'''

network_size=100        #number of  nodes / max ~1000
clump=0.5                #chance of making closest connections
k=int(network_size/2)    #choices of connections out of population
no_connects=3            # approxiate number of starting connections
time=20                 #number of time steps

spec_rate=0.15        #rate of speciation 
ext_rate=0.1            #rate of extinction

mods=False               #node specific rate modifiers?
lineage=False             # lineage specific rate modifiers?
ecosystem_attributes=0  #number of ecosystem attributes
special=False             # special non-clumpy nodes?
special_likely=0      # probability of special nodes
overlap=0        # proportion of occupied neighbour niches to confer competition
runs=10


network='basic network'
version='0.5'


#k and no_connects are related; k (options that connections are chosen from) must be > no_connects
# in order for no_connects to have effect (i.e. cant have 7 connects if only 6 options)
if k>=network_size:
    network_size=k+1
if network_size>k:
    k=network_size-1
if no_connects>k:
    k=no_connects


connections=[]
ecosystem=[random.randint(0,10) for a in range(ecosystem_attributes)]


spec_count=0
ext_count=0

# stop loop if nothig
# identify convergence







'''******************* FUNCTIONS *********************'''

def clumpiness(nodes,clump,network_size,k,no_connects):
    targets=[]
    specials=[]
    for n in nodes:
        connects=[]
        possible_connections = list(range(network_size))
        possible_connections.remove(n)
        options=random.sample(possible_connections,k)
        options.sort()
        s=random.uniform(0,1)
        if s<special_likely:
            specials.append(n)
        while len(connects) < no_connects:
            p=random.uniform(0,1)
            if special==True:
                if s<special_likely:
                    c=random.choice(possible_connections)
                    if c not in connects:
                        connects.append(c)
            elif p<clump:
                    closest=min(options, key=lambda x: abs(x - n))
                    if closest not in connects:
                        connects.append(closest)
                    options.remove(closest)
                    if not options:
                        break  
            else:
                while len(connects) < no_connects and options:
                    c=random.choice(options)
                    if c not in connects:
                        connects.append(c)
                options.remove(c)  
        targets.append(connects)
    return(targets,specials)
        





# SPAWNING THE NETWORK
def init_network(network_size,mods):
    nodes=list(range(network_size))
    targets,specials=clumpiness(nodes,clump,network_size,k,no_connects)
    for no,ta in zip(nodes,targets):
        for t in ta:
            if no not in targets[t]:
                targets[t].append(no)
    sr=[random.random() for _ in range(network_size)]
    er=[random.random() for _ in range(network_size)]
    
    #sp=[random.choice([True,False]) for _ in range(network_size)]
    
    sp=[False for _ in range(network_size)]
    first=random.choice(range(network_size))
    sp[first]=True
    
    counter=0
    tt=[[] for p in sp]
    lssm=[random.random() for _ in range(len(sp))]
    lsem=[random.random() for _ in range(len(sp))]
    for s,(t,v) in enumerate(zip(lssm,lsem)):
        if sp[s]==False:
            lssm[s]=0
            lsem[s]=0
        if sp[s]==True:
            counter=counter+1
            tt[s]=counter
    ln=['' for _ in range(network_size)]
    data={'node id':nodes,'target':targets,'spec. mod':sr,'ext. mod':er,'species':sp,'lineage specific spec. rate mod':lssm,'lineage specific ext. rate mod':lsem ,'lineage':ln,'species tracker':tt}
    df=pd.DataFrame(data)
    df=modifiers(mods,lineage,df)
    dfex=df.explode('target')

    if ecosystem_attributes>0:
        names=['effect','requirement']
        attributes=[random.choice(string.ascii_letters) for _ in range(ecosystem_attributes)]
        attribute_col=[a for a in attributes for _ in range(2)]
        names=(names*((len(attribute_col) // len(names)) + 1))[:len(attribute_col)]       
        starting_effects=[]
        starting_requirements=[]
        stats=[]
        for n in names:
            if n=='effect':
                e=[random.randint(-10,10) for _ in range(network_size)]   
                starting_effects.append(e)
            if n=='requirement':
                mins=[]
                maxs=[]
                for z in range(network_size):
                    r=random.randint(1,10)
                    rmx=random.randint(r,20)
                    mins.append(r)
                    maxs.append(rmx)
                tups=list(zip(mins,maxs))
                starting_requirements.append(tups)
        pairs=list(zip(starting_effects,starting_requirements))
        for p in pairs:
            for i in p:
                stats.append(i)
        stats = list(zip(*stats))
    
        cols=list(zip(attribute_col,names))
        index=pd.MultiIndex.from_tuples(cols, names=['attribute', 'type_']) 
        
        num_columns = len(index)
        ecosystem_engineering_df = pd.DataFrame(stats, columns=index)
        df = pd.concat([df, ecosystem_engineering_df], axis=1)
         
    for i, r in df.iterrows():
        if df.loc[i,'species']==False:
            df.loc[i, df.columns[8:]]=0
    return(df,dfex,counter,specials)
    



def modifiers(mods,lineage,df):
    if mods==False:
        for node, row in df.iterrows():
            df.iloc[node,3]=1
            df.iloc[node,2]=1
    if mods==True:
        pass
    if lineage==False:
        for node, row in df.iterrows():
            df.iloc[node,5]=1
            df.iloc[node,6]=1
            
    if lineage==True:
        pass
    return(df)  
    




def graph(df):
    edges = [(row['node id'], row['target']) for index, row in df.iterrows()]
    G = nx.Graph() #or nx.DiGraph()
    G.add_nodes_from(df['node id'])
    G.add_edges_from(edges)
    
    cmap = plt.get_cmap('Set2')  # Choose a colormap
    norm = mcolors.Normalize(vmin=-2, vmax=4)  # Normalize to 0-1 for boolean values

    color_map = []
    for node in G:
        species_value = df.loc[df['node id'] == node, 'species'].values[0]
        color = cmap(norm(species_value))
        color_map.append(color)

    random_pos = nx.random_layout(G, seed=99)
    pos = nx.spring_layout(G, pos=random_pos)
    nx.draw(G, pos, with_labels=True, node_size=300, node_color=color_map, font_size=10, font_color='black')
    plt.show()
    plt.figure(dpi=200)
    
    # net=Network()
    # net.from_nx(G)
    # net.save_graph("networkx-pyvis.html")
    # HTML(filename="networkx-pyvis.html")
    
#ignore this
# def dynamic_graph(timeline):
#     # PATHPY
#     tn = pp.TemporalNetwork()
#     # PYVIS
#     G = nx.Graph() #or nx.DiGraph()
#     for i,l in enumerate(timeline):
#         edges = [(row['node id'], row['target']) for index, row in l.iterrows()]
        
#         G.add_nodes_from(l['node id'])
#         G.add_edges_from(edges)

#         for e in edges:
#             tn.add_edge(e[0],e[1],i)
        
#     style = {    
#     'ts_per_frame': 1, 
#     'ms_per_frame': 2000,
#     'look_ahead': 2, 
#     'look_behind': 2, 
#     'node_size': 5, 
#     'inactive_edge_width': 2,
#     'active_edge_width': 4, 
#     'label_color' : '#ffffff',
#     'label_size' : '24px',
#     'label_offset': [0,5]
#     }
#     pp.visualisation.export_html(tn, 'C:/Users/Student/OneDrive/Documents/PROJECT/temporal_network.html', **style)
#     #print(tn)

#     net=Network()
#     net.from_nx(G)
#     net.save_graph("networkx-pyvis.html")
#     HTML(filename="networkx-pyvis.html")
    
    
    

def get_neighbours(df,node):            # get list of connections per node
    pairs=[(row['node id'], row['target']) for index, row in df.iterrows()]
    neighbourhoods={}
    for no,ta in pairs:
        if no not in neighbourhoods:
            neighbourhoods[no] = []
        neighbourhoods[no].append(ta)
    neighbourhoods = [[no] + ta for no, ta in neighbourhoods.items()]
    return(neighbourhoods[node])





def neighbour_effect(df,node,neighbours,overlap,ext_rate):          # get proportion of occupied neighbours per node
    small_df = df.groupby('node id').first().reset_index()
    community=[]
    neighbours.remove(node)
    count=(df['node id']==neighbours[0]).sum()
    for n in neighbours:     
        if small_df.iloc[n,4]==True:
            community.append(n)
        if small_df.iloc[n,4]==False:
            pass
        
    # if (len(community)/count)<overlap:
    #     ext_rate=(ext_rate + (ext_rate*(len(community)/count)))
    
    ext_rate=(ext_rate + (ext_rate*(len(community)/count)*(overlap)))
   

    
    return(neighbours,ext_rate)  
   







def speciation(df,dcopy,node,rate,neighbours,counter,spec_count):
    p=random.uniform(0,1)                   # random number
    to_speciate=random.choice(neighbours)   # pick a random neighbour
    #print(dcopy.loc[to_speciate,'species'])
    if p<rate and (dcopy.loc[to_speciate,'species']).all()==False:
        requirement_columns = [col for col in df.columns if isinstance(col, tuple) and col[1] == 'requirement']
        effect_columns=[col for col in df.columns if isinstance(col, tuple) and col[1] == 'effect']
        check=[True]
        for col, e in zip(requirement_columns, ecosystem):
            series=df.loc[to_speciate,[col]]
            req=(series.iloc[0])
            if isinstance(req,tuple):
                if req[0]<=e and req[1]>=e:
                    check.append(True)
                if req[0]>=e and req[1]<=e:
                    check.append(False)
            
        if all(check)==True:                # if all requirements are met then speciate
            dcopy.loc[to_speciate,'species']=True
            spec_count=spec_count+1
            lineage_list=[node,df.loc[node,'lineage']]
            if isinstance(lineage_list, list):
                lineage = ','.join(map(str, lineage_list))
                
            lineage = f"{lineage}"
            dcopy.at[to_speciate, 'lineage'] = lineage
            counter=counter+1
            df.loc[to_speciate,'species tracker']=counter
            dcopy.loc[to_speciate,'lineage specific spec. rate mod']=df.loc[node,'lineage specific spec. rate mod']
            dcopy.loc[to_speciate,'lineage specific ext. rate mod']=df.loc[node,'lineage specific ext. rate mod']
            
            for r,e in zip(requirement_columns,effect_columns):
                dcopy.loc[to_speciate,[r]]=df.loc[node,[r]]
                dcopy.loc[to_speciate,[e]]=df.loc[node,[e]]

            
            
        if any(check)==False:               # if any one requirement not met, dont speciate
            return(dcopy,counter)
            pass
    else:
        pass
    return(dcopy,counter,spec_count)





        

def extinction(dcopy,node,rate,ecosystem,ext_count): 
    extinct=None                # True/False label for each node
    requirement_columns = [col for col in df.columns if isinstance(col, tuple) and col[1] == 'requirement']
    check=[True]
    
    p=random.uniform(0,1)
    if p<rate:
        dcopy.loc[node,'species']=False
        ext_count=ext_count+1
        dcopy.loc[node,'lineage']=''
        dcopy.loc[node,'lineage specific spec. rate mod']=0
        dcopy.loc[node,'lineage specific ext. rate mod']=0
        df.loc[node,'species tracker']=''#
        extinct=True
    if p>rate:
        extinct=False
    
    
    for col, e in zip(requirement_columns, ecosystem):
        series=df.loc[node,[col]]
        req=series.iloc[0]
        if isinstance(req,tuple):
            if req[0]<e and req[1]>e:
                check.append(True)
            if req[0]>e and req[1]<e:
                check.append(False)
  
    if any(check)==False:                # if any requirement not met, go extinct
        dcopy.loc[node,'species']=False
        ext_count=ext_count+1
        dcopy.loc[node,'lineage']=''
        dcopy.loc[node,'lineage specific spec. rate mod']=0
        dcopy.loc[node,'lineage specific ext. rate mod']=0
        df.loc[node,'species tracker']=''
        extinct=True
    
    return(dcopy,extinct,ext_count)
 
 

 
 
 
 
def mods_combine(r,m,l):
    # different options      # adjusted_rate=r*m*l   # r+((m-r)*mod_sensitivity)
    adjusted_rate=r+((m-r)*r)+((l-r)*r)              # currently
    return(adjusted_rate)
 
    
 
    
 
 
    
def ecosystem_engineering(df,ecosystem):
    effect_columns = [col for col in df.columns if col[1] == 'effect']
    effect_sums = df[effect_columns].sum()
    new_ecosystem=[x + y for x, y in zip(ecosystem, effect_sums)]

    return(new_ecosystem)







'''***************** CALLING FUNCTIONS ****************'''

timeline=[]             # initialise list of networks
plt.figure(dpi=200)     # resolution of graph



''' LOOP '''


clumps=[i/10 for i in range(11)]  

folder_name='{} + {}'.format(network,version)

os.makedirs(folder_name, exist_ok=True)

for r in range(10):
    #clump=clumps[r]
    empty=0
    
    df,dfex,counter,specials=init_network(network_size,mods)
    
    file_name='{}{}.csv'.format('clump',clumps[r])
    file_path=os.path.join(folder_name, file_name)
    
    
    settings={'@ SETTINGS':'','Network size':network_size,'Clumpiness':clump,'Special nodes':specials,'Approx. number of connections':no_connects,'Total timesteps':time,'Nodes have specific rate modifiers?':mods,'Species have lineage specific rate modifiers?':lineage,'Number of ecosystem attributes':ecosystem_attributes,'Special nodes allowed?':special,'(If True; likelihood of special nodes?)':special_likely,'Niche overlap':overlap,'Global speciation rate':spec_rate,'Global extinction rate':ext_rate,'Starting ecosystem':ecosystem}
    
    settings_tab=pd.DataFrame(settings.items(), columns=['',''])
    
    settings_tab.to_csv(file_path,mode='w', index=False)

    
    with open(file_path, 'a') as f:
        f.write('\n# Original network\n')
    
    # with open(file_path,'w') as bn:              # write settings and starting network to .csv
    #     bn.write("%s\n" % str('SETTINGS'))
    #     bn.write("%s\n" % str('Network size = {}'.format(network_size))) 
    #     bn.write("%s\n" % str('Clumpiness = {}'.format(clump))) 
    #     bn.write("%s\n" % str('Special nodes = {}'.format(specials))) 
    #     bn.write("%s\n" % str('Approximate number of connections = {}'.format(no_connects)) )
    #     bn.write("%s\n" % str('Total timesteps = {}'.format(time))) 
    #     bn.write("%s\n" % str('Nodes have specific rate modifiers = {}'.format(mods)))
    #     bn.write("%s\n" % str('Species have lineage specific rate modifiers = {}'.format(lineage)))
    #     bn.write("%s\n" % str('Number of ecosystem attributes = {}'.format(ecosystem_attributes))) 
    #     bn.write("%s\n" % str('Special nodes allowed? = {}'.format(special)))
    #     bn.write("%s\n" % str('(If True, likelihood of special nodes?) = {}'.format(special_likely))) 
    #     bn.write("%s\n" % str('Niche overlap = {}'.format(overlap)))
    #     bn.write("%s\n" % str('Global speciation rate = {}'.format(spec_rate)))
    #     bn.write("%s\n" % str('Global extinction rate = {}'.format(spec_rate)))
    #     bn.write("%s\n" % str('')) 
    #     bn.write("%s\n" % str('ORIGINAL NETWORK')) 
    #     bn.write("%s\n" % str('Current ecosystem = {}'.format(ecosystem)))
    
    df.to_csv(file_path,mode='a',index=False)
    graph(dfex)                                     # plot network
    timeline.append(dfex)
    
    for t in range(time):
        spec_count=0
        ext_count=0
        dfex=df.explode('target')
        dcopy=dfex.copy()   
            
        for node,row in df.iterrows():
            node_sm=df.loc[node,'spec. mod']            # identify node specific rate modifiers (1 if False)
            node_em=df.loc[node,'ext. mod']
            lineage_sm=df.iloc[node,5]                  # identify lineage specific rate modifiers (1 if False)
            lineage_em=df.iloc[node,6]
            new_spec_rate=mods_combine(spec_rate,node_sm,lineage_sm)
            new_ext_rate=mods_combine(ext_rate,node_em,lineage_em)
            species_status=df.loc[node,'species']
            neighbours=get_neighbours(dfex,node)
            neighbours,new_ext_rate=neighbour_effect(dfex,node,neighbours,overlap,ext_rate)
            if species_status==True:
                dcopy,counter,spec_count=speciation(df,dcopy,node,new_spec_rate,neighbours,counter,spec_count)
                dcopy,extinct,ext_count=extinction(dcopy,node,new_ext_rate,ecosystem,ext_count)
            
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
       
        graph(dcopy)
        
        ecosystem=ecosystem_engineering(df,ecosystem)   # calculate change in ecosystem to take effect next time step
        t=t+1
        print('TIME = ',t)
        time_title=('TIME = ',t)
        
        
        time_stats={'@ TIME':t,'Current ecosystem':[ecosystem],'Number of speciation events':spec_count,'Number of extinction events':ext_count,'Number of species':species_count}
        time_tab=pd.DataFrame(time_stats.items(),columns=['',''])
        
        time_tab.to_csv(file_path,mode='a', index=False)
       
        with open(file_path, 'a') as f:
            f.write('# Species Network Data\n')
        
        species_network.to_csv(file_path,mode='a', index=False)
        timeline.append(dcopy)                          # add new network to timeline
        if species_count==0:
            empty=empty+1
        if empty==4:
            break

    













