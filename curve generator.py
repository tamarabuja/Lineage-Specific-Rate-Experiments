# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 13:59:03 2024

@author: Student
"""

import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from matplotlib import colormaps
list(colormaps)
import pandas as pd
import os
import ast
import csv
import numpy as np
from scipy import stats
from scipy.stats import pearsonr
from scipy.optimize import curve_fit
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.metrics import r2_score
from scipy.stats import kendalltau
from scipy.stats import t
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression



def exp_model(x, a, b):
    return a *  np.exp(-b * x)

version = 'lineage exp s=e+0.05'
runs = 10

files = []
rc('font', **{'family': 'serif', 'serif': ['Times']})
rc('text', usetex=True)
font = {'fontname': 'serif'}
plt.figure(dpi=1200)


folder_name = 'lineage exp s=e+0.05'

if folder_name == 'lineage exp s=e+0.05':
    graphname = 'Lineage Specific Rate Relationship S=E+0.05'


# graph_file='baseline graphs'


def perlineage(rates, numspecies):
    plrates = []
    rates = rates[1:]
    numspecies = numspecies[1:]
    for i, (x, y) in enumerate(zip(rates, numspecies)):
        if rates[i]!=0 and numspecies[i-1]!=0:
            # print(rates[i],numspecies[i-1])
            plr = rates[i]/numspecies[i-1]
            plrates.append(plr)
        else:
            plrates.append(0)
    # print(plrates)
    return (numspecies, plrates)


def fill(allnums):
    for i in allnums:
        for a in i:
            if len(a) <len(timeline):
                space=len(timeline)-len(a)
                fill=[0]*space
                a+=(fill)
    return(allnums)


def eqm(arr):
    eqmrange = arr[-20:]
    eqmrange = np.array(eqmrange)
    eqm = eqmrange.mean()
    eqm = round(eqm)
    arr = np.array(arr)
    inds = np.where(arr > eqm)[0]
    if inds.size > 0:
        eqmtime = inds[0]
    return (eqm, eqmtime)


for r in range(runs):
    files.append('{}.{}.csv'.format(version, r))


allnums = []
allsr = []
aller = []
alldiv=[]
allpls = [[] for i in range(100)]
allple = [[] for i in range(100)]
allspec = []
allext = []

# file_name = '{}{}.csv'.format(version,run)
files = ['lineage exp s=2e+0.025.1.csv']

indices = []

#directory = 'C:/Users/Student/Downloads/lineage_specific_rates_experiement/lineage exp s=e+0.05'
largerun='C:/Users/Student/Downloads/Large Runs'
parallel='C:/Users/Student/Downloads/lineage exp s=e+0.05.1'
converging ='C:/Users/Student/Downloads/lineage exp s=0.5e+0.075.1(2)'
diverging = 'C:/Users/Student/Downloads/lineage exp s=2e+0.025.2'

run=converging




if run == parallel:
    graph='Parallel'
if run == converging:
    graph='Converging'
if run == diverging:
    graph='Diverging'


#for m, file in enumerate(os.listdir('C:/Users/Student/Downloads/lineage_specific_rates_experiement/lineage exp s=e+0.05')):

for m, file in enumerate(os.listdir('{}'.format(run))):

    print(file)
    file_path = os.path.join('{}'.format(run), file)

# for m,f in enumerate(files):
#     file_path=os.path.join(folder_name, f)
    count = 0

    dataframes = []

    data = {}
    species_net_data={}

    speciations = []
    extinctions = []
    numspecies = []
    timeline = []

    title = None
    with open('{}'.format(file_path), 'r', errors='replace') as file:
        reader = csv.reader(file)
        for line in reader:
            processed = []
            joined_line = ','.join(line).strip()
            if not line:
                continue

            if line[0].startswith('#'):
                joined_line = ','.join(line).strip()
                title = joined_line[1:].strip()
                # print(title)
                if title not in data:
                    data[title] = []
                continue

            if line[0].startswith('@'):
                joined_line = ','.join(line).strip()
                title = joined_line[1:].strip()
                # print(title)
                if title not in data:
                    data[title] = []
                continue

            if title:
                for v in line:
                    v = v.strip()
                    if v.startswith('[') and v.endswith(']'):
                        v = ast.literal_eval(v)

                    processed.append(v)
                if processed:
                    data[title].append(processed)

    for i, (key, value) in enumerate(data.items()):
        if 'TIME' in key:
            # print(key)
            time_integer = int(key.split(',')[1])
            # print(time_integer)
            timeline.append(time_integer)
            speciations.append(value[0][1])
            extinctions.append(value[1][1])
            numspecies.append(value[2][1])
    count = count+1

    #TOTAL
    speciations = [int(s) for s in speciations]
    extinctions = [int(s) for s in extinctions]
    numspecies = [int(s) for s in numspecies]
    # print(numspecies)
    # print(speciations)
    # print(extinctions)




    # COUNT EACH LINEAGE IN EACH TIMESTEP
    #print(taxa)
    iterdata=(data['Species Network Data'])
    #print(iterdata)
    table=[]
    subtables=[]
    for r in (iterdata):
        if 'node id' in r:
            if table:
                subtables.append(table)
                table=[]
        if r and r[0]!='':
            table.append(r)
    if table:
        subtables.append(table)
    
    svalue_counts={}   
    excounts={}
    #print(subtables[0])
    o=0
    for st in subtables:
        o=o+1
        specs=[sub[3] for sub in st]
        extes=[sub[4] for sub in st]
        #print(specs)
        spec_count=pd.Series(specs[1:]).value_counts(sort=False)
        ext_counts=pd.Series(extes[1:]).value_counts(sort=False)
        # print()
        # print(spec_count)
        for value, count in spec_count.items():
            if value in svalue_counts:
                svalue_counts[value].append(count)
            else:
                svalue_counts[value] = [count]      
        for x, counte in ext_counts.items():
            if x in excounts:
                excounts[x].append(counte)
            else:
                excounts[x] = [counte]      
       
    #print(value_counts)
    exr=list(excounts.keys())
    taxa=list(svalue_counts.keys())
    nums=list(svalue_counts.values())
    extinct=[[] for _ in nums]
  
    for n,m  in zip(nums,extinct):
        for i in range(1, len(n)):
            if n[i] < n[i - 1]: 
                ex = n[i - 1] - n[i]  
                m.append(ex)  
            elif n[i]>=n[i-1]:
                m.append(0)  
   
    
    
    # for y,t in zip(nums,extinct):
    #     div,pler=perlineage(y,t)
        
    
    
    allnums.append(nums)
    allext.append(extinctions)
    allspec.append(numspecies)

    allnums=fill(allnums)
    #allext=fill(allext)
               
    
        #plt.plot(timeline,nums[i], label=taxa[i])
    
   
#allnums=[arr for arr in allnums if len(arr) == 6 and all(len(sub_arr) == 1001 for sub_arr in arr)]


allnums = [np.array(item) for item in allnums]
allext = [np.array(item) for item in allext]
allspec = [np.array(item) for item in allspec]


stackext=np.stack(allext)
stacknum = np.stack(allnums)
stackspec = np.stack(allspec)

meannum = np.mean(stacknum, axis=0)
meanext=np.mean(stackext, axis=0)
meansp=np.mean(stackspec, axis=0)

stdnum = np.std(stacknum, axis=0)
semnum = stdnum / np.sqrt(stacknum.shape[0])

divs=[]
plers=[]
div,pler=perlineage(extinctions,numspecies)
divs.append(div)
plers.append(pler)
    

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']


for f,y in enumerate(meannum):
    # r_value, p_value = pearsonr(timeline, y)
    # r_squared = r_value ** 2  # Use r_squared for clarity
 
    # if 0.1 <= p_value <= 1:
    #     # For p-values within this range, use standard decimal notation
    #     p_label = rf"${p_value:.2f}$"
    # else:
    #     # Otherwise, use scientific notation
    #     p_value_str = "{:.2e}".format(p_value)
    #     p_base, p_exp = p_value_str.split('e')
    #     p_exp = int(p_exp)  # Convert exponent to integer

    #     # Format p_label to exclude 10^0 if exponent is zero
    #     if p_exp == 0:
    #         p_label = rf"${float(p_base):.2f}$"
    #     else:
    #         p_label = rf"${p_base} \times 10^{{{p_exp}}}$"
    plt.plot(timeline,y,color=CB_color_cycle[f],label=rf"S={float(taxa[f]):.4f}, E={float(exr[f]):.4f}")
    plt.fill_between(timeline, meannum[f] - semnum[f], meannum[f] + semnum[f], alpha=0.15,color=CB_color_cycle[f])
    
 
plt.title('Rate of Speciation, {} Rates'.format(graph))
plt.xlabel('Time')
plt.ylabel('Mean No. of Species')
plt.legend(fontsize='small')
plt.show()
plt.figure(dpi=1200)

# plt.stackplot(timeline,*meannum,labels=[f"{float(t):.4f}" for t in taxa])
# plt.xlim(0,500)
# #plt.yscale("log")   
# plt.title('Rate of Speciation, {} Rates'.format(graph))
# plt.xlabel('Time')
# plt.ylabel('Mean No. of Species')
# plt.legend(loc='upper left')
# plt.show()
# plt.figure(dpi=1200)

timeline=timeline[1:]

start=50
timeline=timeline[start:]
pler=pler[start:]


stdext = np.std(stackext, axis=0)
semext = stdext / np.sqrt(stackext.shape[0])

plt.plot(timeline,pler,linewidth=0.3,color='#ff7f00') #marker='.',s=15)
#plt.fill_between(timeline[1:],pler - stdext[1:], pler + stdext[1:], alpha=0.15)

plt.title('Per-Lineage Extinction Probability, {} Rates'.format(graph))
plt.xlabel('Time')
plt.ylabel('Extinction Probability')


'''LINEAR REGRESSION '''
e_slope, e_intercept, e_r_value, e_p_value, e_std_err = stats.linregress(timeline, pler)
x_range = np.linspace(min(timeline), max(timeline), 1000)
e_y_pred = e_slope * x_range + e_intercept
esmoothed = lowess(pler, timeline, frac=0.2)
ex_lowess = esmoothed[:, 0]
ey_lowess = esmoothed[:, 1]
e_r_value = r2_score(pler, np.interp(timeline, ex_lowess, ey_lowess))
e_label_text = (r'P = {:1e}, $R^2$ = {:.4f}'.format(e_p_value, e_r_value**2))
R2=e_r_value**2

def split(e_p_value):
    coefficient = f"{ e_p_value:.6g}".split('e')[0]  
    exponent = int(f"{ e_p_value:.6e}".split('e')[1])  
    return(coefficient,exponent)

coefficient, exponent=split(e_p_value)
# Create the label using MathTeX formatting
label = fr"$P = {coefficient} \times 10^{{{exponent}}},\ R^2 = {R2:.4f}$"
#plt.plot(x_range, e_y_pred, color='#e41a1c', linestyle='--', linewidth=2,label=label)


''' POLYNOMIAL REGRESSION '''
x=timeline
y=pler
x=np.array(x)
x = x.reshape(-1, 1)


poly = PolynomialFeatures(degree=1)  # Experiment with different degrees
x_poly = poly.fit_transform(x)
# Fit the polynomial regression model
model = LinearRegression()
model.fit(x_poly, y)

# Predict y values for the fitted line
y_pred = model.predict(x_poly)

r_squared = r2_score(y, y_pred)
n = len(y)  # Number of observations
p = x_poly.shape[1] - 1  # Number of predictors (exclude intercept)
adjusted_r_squared = 1 - ((1 - r_squared) * (n - 1) / (n - p - 1))

slope_p_value = kendalltau(x, y).pvalue  # Mann-Kendall test for trend

coefficient,exponent=split(slope_p_value)
label = fr"$P = {coefficient} \times 10^{{{exponent}}},\ $Adjusted $R^2 = {adjusted_r_squared:.4f}$"
plt.plot(x, y_pred, color='red', label=label)



def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')

def r_squared(original, smoothed):
    mask = ~np.isnan(smoothed)  # Exclude NaN values from calculation
    original_valid = original[mask]
    smoothed_valid = smoothed[mask]
    ss_res = np.sum((original[mask] - smoothed[mask]) ** 2)  # Residual sum of squares
    ss_tot = np.sum((original[mask] - np.mean(original[mask])) ** 2)  # Total sum of squares
    return 1 - (ss_res / ss_tot)

window_size = 50
pler = np.array(pler)  # Convert original data to NumPy array

y_moving_avg = moving_average(pler, window_size)
y_moving_avg_full = np.concatenate((np.full(window_size - 1, np.nan), y_moving_avg))
y_moving_avg_full = np.array(y_moving_avg_full)  # Ensure smoothed data is a NumPy array

r2_moving_avg = r_squared(pler, y_moving_avg_full)

#plt.plot(timeline[1:][window_size - 1:], y_moving_avg, label=f"Moving Average ($R^2 = {r2_moving_avg:.3f}$)", linewidth=1, color="#e41a1c", linestyle='--')

plt.xlim(50,1000)

#plt.yscale('log')
#plt.xscale('log')
plt.legend()
plt.show()
plt.figure(dpi=1200)


# standard error envelope
# stakced bars
# use lynter to clean code







