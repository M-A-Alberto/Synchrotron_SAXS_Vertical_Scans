# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 18:04:12 2024

@author: amart
"""



import numpy as np


import pandas as pd

import seaborn as sns

import matplotlib.pyplot as plt

import os 


from sklearn.metrics import auc

from scipy.optimize import curve_fit

sns.set_theme()
sns.set_style("ticks")
sns.set_context("paper",font_scale=1)
sns.set_palette("bright")

import warnings
warnings.filterwarnings("error",category = RuntimeWarning)

def recta(x,m,a):
    
    y = 10**a*x**m
    
    return y


os.chdir(r"XXXXXXXXXXX")


scans = [i for i in os.listdir("Processed") if "snapascan" in i]

blanks = [i for i in os.listdir("./Processed/LaVue1D") if ("pilatus" in i) and ("Blank" in i)]

bkg_channel = pd.read_table(r"Processed/LaVue1D/pilatus_NPs_empty_000_0000.dat",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])

stz_center = -27.534
sx_center = -0.184

table_areas = pd.DataFrame(columns=["Deltaz"])


def bkgs(scans):
    os.makedirs("./z scans both/SAXS/Time evolution/Before_bkg_subtraction",exist_ok=True)
    os.makedirs("./z scans both/SAXS/Time evolution/After_bkg_subtraction",exist_ok=True)
    
    
    t = 0
    for n,scan in enumerate(scans):
    
        positions = pd.read_table(f"raw/{scan}/{scan}.dat",skiprows = 6,header= None,sep=" ",skipfooter=1,engine="python", names = ["Pt_No",  "stz",  "uxtimer",  "pilatus_timer",  "rayonix_timer",  "pilatus_roi",  "rayonix_roi",  "tfg_timer",  "tfg_ch1",  "tfg_ch2",  "tfg_ch3",  "tfg_ch4", "current",  "dt"])
    
        files = [i for i in os.listdir(f"Processed/{scan}/data/LaVue1D/") if "pilatus" in i]
        
        bkg_matrigel = blanks[n]
        
        bkg_matrigel = pd.read_table(f"./Processed/LaVue1D/{bkg_matrigel}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
        
        
        colors = plt.cm.gnuplot(np.linspace(0,1,len(files)))
        
        print(positions)
        
        fig, (ax1,ax2) = plt.subplots(1,2)
        k=0
        for file in files:
            
            No = file[:-4].split("_")[-1]
            
            data = pd.read_table(f"Processed/{scan}/data/LaVue1D/{file}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
            
            stz = positions.loc[positions.loc[:,"Pt_No"]==int(No),"stz"]
            deltaz = np.abs(float(stz.iloc[0]))-np.abs(positions.loc[0,"stz"])
            deltaz = abs(round(deltaz*1000))
            
            if deltaz in [61,63]:
                deltaz = 62
            
            elif deltaz == 124:
                deltaz = 125
                
            elif deltaz == 187:
                deltaz = 188
            
            elif deltaz == 249:
                deltaz = 250
            
            elif deltaz == 312:
                deltaz = 313
            
            elif deltaz == 374:
                deltaz = 375
                
            elif deltaz == 437:
                    deltaz = 438
            
            elif deltaz in [499,501]:
                deltaz = 500
            if int(No)<=1:
                
                ax1.plot(data.q,data.I,color=colors[k],label=str(deltaz))
                
                #data.I = data.I-bkg_channel.I
                
                
            else:
                ax2.plot(data.q,data.I,color=colors[k],label=str(deltaz))
                #data.I = data.I-bkg_matrigel.I
            k+=1
            
            
            
        ax1.plot(bkg_channel.q,bkg_channel.I,"g",label="Channel")
            
        ax2.plot(bkg_matrigel.q,bkg_matrigel.I,"g",label="Matrigel")
        
        ax1.legend()
        ax2.legend()
        
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xscale("log")
        
        ax1.set_title("Channel")
        ax2.set_title("Matrigel")
        
        ax1.set_ylabel("Integrated Scattered Intensity (a.u.)")
        fig.text(0.5, 0.04, 'q (nm$^{-1}$)', ha='center')
    
        
        plt.suptitle(f"Time = {t} mins")
        
        plt.tight_layout()
        
        plt.savefig(f"./z scans both/SAXS/Time evolution/Before_bkg_subtraction/Time_{t}_mins.tif",dpi=300,bbox_inches="tight")
        
        plt.show()
        plt.close()
        
        t+=5
       
    t = 0
    
       
    for n,scan in enumerate(scans):
    
        positions = pd.read_table(f"raw/{scan}/{scan}.dat",skiprows = 6,header= None,sep=" ",skipfooter=1,engine="python", names = ["Pt_No",  "stz",  "uxtimer",  "pilatus_timer",  "rayonix_timer",  "pilatus_roi",  "rayonix_roi",  "tfg_timer",  "tfg_ch1",  "tfg_ch2",  "tfg_ch3",  "tfg_ch4", "current",  "dt"])
    
        files = [i for i in os.listdir(f"Processed/{scan}/data/LaVue1D/") if "pilatus" in i]
        
        bkg_matrigel = blanks[n]
        
        bkg_matrigel = pd.read_table(f"./Processed/LaVue1D/{bkg_matrigel}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
        
        
        colors = plt.cm.gnuplot(np.linspace(0,1,len(files)))
        
        
        print(positions)
        
        fig, (ax1,ax2) = plt.subplots(1,2)
        k=0
        for file in files:
            
            No = file[:-4].split("_")[-1]
            
            data = pd.read_table(f"Processed/{scan}/data/LaVue1D/{file}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
            
            stz = positions.loc[positions.loc[:,"Pt_No"]==int(No),"stz"]
            deltaz = np.abs(float(stz.iloc[0]))-np.abs(positions.loc[0,"stz"])
            deltaz = abs(round(deltaz*1000))
            
            if deltaz in [61,63]:
                deltaz = 62
            
            elif deltaz == 124:
                deltaz = 125
                
            elif deltaz == 187:
                deltaz = 188
            
            elif deltaz == 249:
                deltaz = 250
            
            elif deltaz == 312:
                deltaz = 313
            
            elif deltaz == 374:
                deltaz = 375
                
            elif deltaz == 437:
                    deltaz = 438
            
            elif deltaz in [499,501]:
                deltaz = 500
            
            
            os.makedirs(f"z scans both/SAXS/Data without bkg/Time_{t}_mins",exist_ok=True)
            
            if int(No)<=1:
                
                data.I = data.I-bkg_channel.I
                ax1.plot(data.q,data.I,color=colors[k],label=str(deltaz))
                
                #
                
                
            else:
                data.I = data.I-bkg_matrigel.I
                ax2.plot(data.q,data.I,color=colors[k],label=str(deltaz))
                #
            
            
            data.to_csv(f"z scans both/SAXS/Data without bkg/Time_{t}_mins/{deltaz}_um.txt",index=False,sep="\t")
            
            k+=1
            
            
            
        
        ax1.legend()
        ax2.legend()
        
        ax1.set_yscale("log")
        ax1.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xscale("log")
        
        ax1.set_title("Channel")
        ax2.set_title("Matrigel")
        
        ax1.set_ylabel("Integrated Scattered Intensity (a.u.)")
        fig.text(0.5, 0.04, 'q (nm$^{-1}$)', ha='center')
    
        
        plt.suptitle(f"Time = {t} mins")
        
        plt.tight_layout()
        
        plt.savefig(f"./z scans both/SAXS/Time evolution/After_bkg_subtraction/Time_{t}_mins.tif",dpi=300,bbox_inches="tight")
        
        plt.show()
        plt.close()
        
        t+=5
    
    
    return

def gauss(x,A,mu,s):
    
    y = A*np.exp(-(x-mu)**2/(2*s**2))
    
    return y 
#bkgs(scans)

t = 0

table_areas = pd.DataFrame(columns=["Deltaz"])

i=0

bkg_channel.I = (bkg_channel.I)/bkg_channel.I.iloc[150]



os.makedirs(r"./z scans both/SAXS/Concentration evolution",exist_ok=True)
for n,scan in enumerate(scans):

    positions = pd.read_table(f"raw/{scan}/{scan}.dat",skiprows = 6,header= None,sep=" ",skipfooter=1,engine="python", names = ["Pt_No",  "stz",  "uxtimer",  "pilatus_timer",  "rayonix_timer",  "pilatus_roi",  "rayonix_roi",  "tfg_timer",  "tfg_ch1",  "tfg_ch2",  "tfg_ch3",  "tfg_ch4", "current",  "dt"])

    files = [i for i in os.listdir(f"Processed/{scan}/data/LaVue1D/") if "pilatus" in i]
    
    bkg_matrigel = blanks[n]
    
    bkg_matrigel = pd.read_table(f"./Processed/LaVue1D/{bkg_matrigel}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
    
    bkg_matrigel.I = (bkg_matrigel.I)/bkg_matrigel.I.iloc[150]
    colors = plt.cm.gnuplot(np.linspace(0,1,len(files)))
    
    os.makedirs(f"./z scans both/SAXS/Time evolution/Time {t} mins",exist_ok=True)
    
    areas = []
    
    z = []
    
    
    
    if len(files) != 9:
        print(f"{scan} requires attention")
        
        
    for file in files:
        
        fig, (ax1,ax2) = plt.subplots(1,2)
        No = file[:-4].split("_")[-1]
        
        data = pd.read_table(f"Processed/{scan}/data/LaVue1D/{file}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
        
        stz = positions.loc[positions.loc[:,"Pt_No"]==int(No),"stz"]
        deltaz = np.abs(float(stz.iloc[0]))-np.abs(positions.loc[0,"stz"])
        deltaz = abs(round(deltaz*1000))
        
        
        
        
        if deltaz in [61,63]:
            deltaz = 62
        
        elif deltaz == 124:
            deltaz = 125
            
        elif deltaz == 187:
            deltaz = 188
        
        elif deltaz == 249:
            deltaz = 250
        
        elif deltaz == 312:
            deltaz = 313
        
        elif deltaz == 374:
            deltaz = 375
            
        elif deltaz == 437:
                deltaz = 438
        
        elif deltaz in [499,501]:
            deltaz = 500
            
            
            
        os.makedirs(f"./z scans both/SAXS/Time evolution/Deltaz {deltaz} um",exist_ok=True)
        
        
        data.I = (data.I)/data.I.iloc[150]#-np.mean(data.I.iloc[-100:])
        
        
        
        if int(No)<=1:
            
            data.I = data.I#-bkg_channel.I
            #ax1.plot(bkg_channel.q,bkg_channel.I,label="Channel",color="g")
                  
            
        else:
            data.I = data.I#-bkg_matrigel.I
            #ax1.plot(bkg_matrigel.q,bkg_matrigel.I,label="Matrigel",color="g")
            
        #ROI_data = data.loc[(data.q >= 0.15) & (data.q<=0.43),:].reset_index()
        ROI_data = data.loc[(data.q >= 0.43) & (data.q<=0.7) & (data.I>0),:].reset_index()
        #ROI_data.I = ROI_data.I/ROI_data.I.iloc[0]
        ax1.plot(data.q,data.I,label="Data",color="k")
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        
        #Oax1.set_ylim([0.1,200])
        ax1.set_xlabel("q (nm$^{-1}$)")
        ax1.set_ylabel("Integrated Scattered Intensity (a.u.)")
        try:
            
            y1 = np.log10(ROI_data.iloc[0,2])
            y2 = np.log10(ROI_data.iloc[-1,2])
            
            x1 = np.log10(ROI_data.iloc[0,1])
            x2 = np.log10(ROI_data.iloc[-1,1])
            
            m = (y2-y1)/(x2-x1)
            
            a = y2-m*x2
                
            
            area = auc(ROI_data.q,ROI_data.I-recta(ROI_data.q,m,a))
            
            
            
            ax1.plot(ROI_data.q,recta(ROI_data.q,m,a),"g",label="Background")
            
            
                
           
            ax1.fill_between(ROI_data.q,ROI_data.I,recta(ROI_data.q,m,a),color = "r",alpha=0.3)
            
            #plt.ylim([5e-6,1.5])
                
            ax2.plot(ROI_data.q,ROI_data.I-recta(ROI_data.q,m,a))
            y_fit = np.array(ROI_data.I-recta(ROI_data.q,m,a))
            x_fit = np.array(ROI_data.q)
        
            pars, cov = curve_fit(gauss,x_fit,y_fit,p0=[0.1,0.5,0.1])
            
            ax2.plot(ROI_data.q,gauss(ROI_data.q,*pars))
            
            area = pars[0]
        
        except:
           
            print(f"File {file} failed")
            area = 0
        
        ax2.set_xlabel("q (nm$^{-1}$)")
        ax1.legend(loc="upper right")
        plt.suptitle(f"$\Delta$z = {deltaz} $\mu$m, Time = {t} mins")
        #plt.savefig(f"./z scans both/SAXS/Time evolution/Time {t} mins/{deltaz}_um.tif",dpi=100,bbox_inches="tight")
        #plt.savefig(f"./z scans both/SAXS/Time evolution/Deltaz {deltaz} um/Time {t} mins.tif",dpi=100,bbox_inches="tight") 
        
        plt.show()
        plt.close()
        
        z.append(deltaz)
        areas.append(area)
    
        table_areas.loc[i,"Deltaz"] = deltaz
        table_areas.loc[i,"Areas"] = area
        table_areas.loc[i,"Time (mins)"] = t
        
        
        i+=1
            
    areas = np.array(areas)
    
    areas= areas/areas[0]
    
    plt.plot(z,areas,".-")
    
    ax = plt.gca()
    
    ylims = ax.get_ylim()
    
    ax.fill_betweenx(ylims,-50,100,hatch= "O",alpha=0.3,color="peru")
    ax.fill_betweenx(ylims,100,550,hatch="/",alpha=0.3,color="pink")
    
    ax.set_ylim(ylims)
    
    plt.xlabel("Distance ($\mu$m)")
    plt.ylabel("Normalized Concentration")
    plt.title(f"Time = {t} mins")
    
    #plt.savefig(f"./z scans both/SAXS/Concentration evolution/Time {t} mins.tif",dpi=100,bbox_inches="tight")
    
    plt.show()
    plt.close()
    
    
    
    
    t+=5
    



for t in set(table_areas.loc[:,"Time (mins)"]):
    
    
    select = table_areas.loc[table_areas.loc[:,"Time (mins)"]==t,:]

    table_areas.loc[table_areas.loc[:,"Time (mins)"]==t,"Areas"]= table_areas.loc[table_areas.loc[:,"Time (mins)"]==t,"Areas"]/float(select.loc[select.loc[:,"Deltaz"]==0,"Areas"].iloc[0])
    
    print(table_areas.loc[table_areas.loc[:,"Time (mins)"]==t,:])



sns.lineplot(data=table_areas,x="Deltaz",y = "Areas",hue="Time (mins)",palette="bright",marker="o")
ax = plt.gca()
ylims = ax.get_ylim()
xlims = ax.get_xlim()

ax.fill_betweenx(ylims,-50,100,hatch= "O",alpha=0.3,color="peru")
ax.fill_betweenx(ylims,100,550,hatch="/",alpha=0.3,color="pink")

ax.set_ylim(ylims)
ax.set_xlim(xlims)
plt.xlabel("Distance ($\mu$m)") 
plt.ylabel("SAXS Area (a. u.)")      
plt.legend(ncol=3)

 
plt.savefig("./z scans both/SAXS/Time Evolution/Distance.tif",dpi=300,bbox_inches="tight")        
        
plt.show()
plt.close()



sns.lineplot(data=table_areas,x="Time (mins)",y = "Areas",hue="Deltaz",palette="bright",marker="o")



plt.ylabel("SAXS Area (a. u.)")       
plt.savefig("./z scans both/SAXS/Time Evolution/Time.tif",dpi=300,bbox_inches="tight")        
        
plt.show()
plt.close()







blank_area = []
time = []

t = 0

for blank in blanks:
    
    data = pd.read_table(f"./Processed/LaVue1D/{blank}",skiprows = 23,header=None,sep = "\s+",names = ["q","I","sigma"])
    
    plt.plot(data.q,data.I,label=f"{t} mins")
     
    ROI_data = data.loc[(data.q >= 0.43) & (data.q<=0.7),:].reset_index()
    
    y1 = np.log10(ROI_data.iloc[0,2])
    y2 = np.log10(ROI_data.iloc[-1,2])
    
    x1 = np.log10(ROI_data.iloc[0,1])
    x2 = np.log10(ROI_data.iloc[-1,1])
    
    m = (y2-y1)/(x2-x1)
    
    a = y2-m*x2
        
    
    area = auc(ROI_data.q,ROI_data.I-recta(ROI_data.q,m,a))
    
    
    blank_area.append(area)
    time.append(t)
    
    
    t+=5


plt.xscale("log")
plt.yscale("log")

plt.legend(ncol=3)
plt.xlabel("q (nm$^{-1}$)")
plt.ylabel("Integrated Scattered Intensity (a.u.)")
plt.savefig("./z scans both/SAXS/Time Evolution/Blanks.tif",dpi=300,bbox_inches="tight")
plt.show()
plt.close()


blank_areas = pd.DataFrame()

blank_areas.loc[:,"Time (mins)"] = time
blank_areas.loc[:,"Blank Area (a. u.)"] = blank_area

sns.lineplot(data = blank_areas,x = "Time (mins)",y = "Blank Area (a. u.)",palette = "tab10",marker="o")

plt.savefig("./z scans both/SAXS/Time Evolution/Blanks_areas.tif",dpi=300,bbox_inches="tight")





