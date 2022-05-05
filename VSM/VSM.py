#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VSM - Version 1.0

   Volcanic and Seismic source Modelling

   Author:  Elisa Trasatti, elisa.trasatti@ingv.it
   Istituto Nazionale di Geofisica e Vulcanologia - Rome (Italy)

   Thanks to: M. Sambridge (NA inversion)
              N. Piana Agostinetti (first fortran version)
              And others for python codes (from GitHub and PyPI)

   Created:      2007 
   Last update:  May 2022

   License:  E. Trasatti, covered by GNU-GPL License

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

      * Redistributions of source code must retain the above copyright 
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright 
 	    notice, this list of conditions and the following disclaimer in 
 	    the documentation and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPY RIGHT OWNER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
"""

import os

#import sys
#sys.path.append('path-to-VSM-folder')
import VSM_search as NAsearch
import VSM_forward as forward

import numpy as np
np.random.seed(42)
import pandas as pd
import shapefile
import datetime

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox

from matplotlib import pyplot as plt
from matplotlib import rcParams
import corner
import emcee

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
plt.rcParams["font.family"] = "Calibri"
plt.rcParams.update({'font.size':14})

# -----------------------------------------------------------------------------
# INPUT DATA

ff_sar = []
ff_sar.append('None')
ff_gps = 'None'
ff_lev = 'None'
ff_edm = 'None'
ff_tlt = 'None'
ff_srn = 'None'

n_sources = 0
mysources = []
okada_mode = np.full(10,None)
bounds = []
labels = []
parflag = []
bounds_true = []
labels_true = []

# GUI ONLY ---------------
def read_inout_folder():
    global fold_inout
    fold_inout = filedialog.askdirectory(title = "Select Folder")

def final_check():
    if 'mu' not in globals():
        read_elast(5.e9,0.25)
    if 'w_sar' not in globals():
        read_weights(1.,1.,1.,1.,1.)
    if 'fold_inout' not in globals():
        window = tk.Tk()
        window.wm_withdraw()
        window.geometry("1x1+200+200")
        txt = 'Forgot to choose In/Out Folder'
        messagebox.showerror(title="Error",message=txt,parent=window)
        read_inout_folder()
        window.destroy()

    ff_in = os.path.join(fold_inout,'VSM_input.txt')
    fguiin = open(ff_in, "w")

    fguiin.write(fold_inout+'\n')
    
    print('GUI',len(ff_sar),ff_sar)    
    try:
        if len(ff_sar)<2:
            fguiin.write(ff_sar[0]+'\n')
        else:
            ff_sar_temp = ''
            for i in range(len(ff_sar)):
                ff_sar_temp += ff_sar[i]+' '
            fguiin.write(ff_sar_temp+'\n')
    except:
        fguiin.write(ff_sar[0]+'\n')
        
    fguiin.write(ff_gps+'\n'+ff_lev+'\n'+ff_edm+'\n'+ff_tlt+'\n'+ff_srn+'\n')
    fguiin.write(str(w_sar)+'\n'+str(w_gps)+'\n'+str(w_lev)+'\n'+str(w_edm)+'\n'+str(w_tlt)+'\n'+str(w_srn)+'\n')
    fguiin.write(('%.2e\n' % mu)+str(ni)+'\n')
    fguiin.write(str(len(mysources))+'\n')
    kk = 0
    for i_sorg in range(len(mysources)):
        if(mysources[i_sorg]['index']!= 5):
            fguiin.write(str(mysources[i_sorg]['index'])+'\n')
        else:
            fguiin.write(str(mysources[i_sorg]['index'])+' '+okada_mode[i_sorg]+'\n')
           
        for k in range(mysources[i_sorg]['param_n']):
            row = str(list(bounds[kk])[0])+'\t'+str(list(bounds[kk])[1])+'\t'+mysources[i_sorg]['param_l'][k]+'\n'
            kk +=1            
            fguiin.write(row)
    fguiin.write(str(ch_opt)+'\n')
    fguiin.write(str(samp1)+' '+str(samp2)+'\n')
    fguiin.write(str(samp3)+'\n')
    fguiin.write(str(nskip)+'\n')
    fguiin.close()



def read_sar_file():
    for i in range(10):
        f = filedialog.askopenfilename(title = "Select InSAR File",filetypes = (("ascii files","*.txt"),("csv files","*.csv"),("shapefiles","*.shp")))
        if(f!=''): 
            if(i == 0):
                ff_sar[0] = f
            else:
                ff_sar.append(f)
        else:
            break
        window = tk.Tk()
        window.wm_withdraw()
        window.geometry("1x1+200+200")
        title='VSM - SAR Data Input File'
        txt = 'Select another SAR file?'
        MsgBox = messagebox.askquestion (title,txt,icon = 'warning')
        window.destroy()
        if MsgBox == 'no':
            break

def read_gps_file():
    global ff_gps
    ff_gps = filedialog.askopenfilename(title = "Select GPS File",filetypes = (("ascii files","*.txt"),("csv files","*.csv"),("shapefiles","*.shp")))
def read_lev_file():
    global ff_lev
    ff_lev = filedialog.askopenfilename(title = "Select LEV File",filetypes = (("ascii files","*.txt"),("csv files","*.csv"),("shapefiles","*.shp")))

def read_edm_file():
    global ff_edm
    ff_edm = filedialog.askopenfilename(title = "Select EDM File",filetypes = (("ascii files","*.txt"),("csv files","*.csv"),("shapefiles","*.shp")))

def read_tlt_file():
    global ff_tlt
    ff_tlt = filedialog.askopenfilename(title = "Select TILT File",filetypes = (("ascii files","*.txt"),("csv files","*.csv"),("shapefiles","*.shp")))

def read_srn_file():
    global ff_srn
    ff_srn = filedialog.askopenfilename(title = "Select STRAIN File",filetypes = (("ascii files","*.txt"),("csv files","*.csv"),("shapefiles","*.shp")))


# ALL ----------------
def read_elast(mm,nn):
    global mu,ni
    mu = float(mm)
    ni = float(nn)

def read_Okada_mode(a):
    okada_mode[n_sources] = a

def read_optim(ch,sa,it,re):
    global ch_opt, samp1, samp2,samp3
    ch_opt = int(ch)
    if(ch_opt != 0 and ch_opt != 1):
            print('\nError! Inversion tool can be NA (choice 0) or BI (choice 1)')
            print('Setting the inversion tool to NA')
            ch_opt = 0
    
    samp1  = int(sa)
    samp3  = int(it)
    if(re != None):
        samp2  = int(re)
    else:
        if(ch_opt==0): samp2 = int(samp1/2)
        if(ch_opt==1): samp2 = 10 # bin dimension

def read_params(m1,m2,opt):
    global labels, bounds,parflag
    global bounds_true, labels_true
    global n_sources
    n_sources += 1
    
    mysources.append(source_info(opt))
    ll = mysources[n_sources-1]['param_l']
            
    if(n_sources == 1): labels += ll
    else: labels+= [l+'_'+str(n_sources) for l in ll]
    
    for i_param in range(mysources[n_sources-1]['param_n']):
        min_val = float(m1[i_param])
        max_val = float(m2[i_param])
        if(labels[i_param] == 'depth' and (max_val < 0. or min_val < 0)):
            print('\nError! Depths are considered positive!')
            print('Considering depth positive for you!')
            if(max_val < 0.):
                max_val = -max_val
            if(min_val < 0.):
                min_val = -min_val
        if(max_val < min_val):
            print('\nError! Max value smaller than min value for parameter:',labels[i_param])
            print('Switching max and min values for you!')
            app     = max_val
            max_val = min_val
            min_val = app
        
        bb = (min_val,max_val)
        bounds.append(bb)
       
        if(min_val != max_val):
            parflag.append(True)
            bounds_true.append(bounds[-1])
            labels_true.append(labels[len(bounds)-1])
        else:
            parflag.append(False)

def read_plot(brn_in):
    global plotyn,nskip
    plotyn= True
    nskip = int(brn_in)
    if(nskip == -1): plotyn = False
        
def read_shapefile(shp_path):
	sf = shapefile.Reader(shp_path)
	fields = [x[0] for x in sf.fields][1:]
	records = sf.records()
	df = pd.DataFrame(columns=fields, data=records)
	return df

def read_weights(s,g,l,e,t,r):
    global w_sar,w_gps,w_lev,w_edm,w_tlt,w_srn
    w_sar = float(s)
    w_gps = float(g)
    w_lev = float(l)
    w_edm = float(e)
    w_tlt = float(t)
    w_srn = float(r)
    
    if(ff_sar[0] == 'None' or ff_sar[0] == ''): w_sar = 0.
    if(ff_gps == 'None'): w_gps = 0.
    if(ff_lev == 'None'): w_lev = 0.
    if(ff_edm == 'None'): w_edm = 0.
    if(ff_tlt == 'None'): w_tlt = 0.
    if(ff_srn == 'None'): w_srn = 0.
    
    if((ff_sar[0] != 'None' and ff_sar[0] != '') and w_sar == 0.): w_sar = 1.
    if(ff_gps != 'None' and w_gps == 0.): w_gps = 1.
    if(ff_lev != 'None' and w_lev == 0.): w_lev = 1.
    if(ff_edm != 'None' and w_edm == 0.): w_edm = 1.
    if(ff_tlt != 'None' and w_tlt == 0.): w_tlt = 1.
    if(ff_srn != 'None' and w_srn == 0.): w_srn = 1.
    
    w_tot = w_sar + w_gps + w_lev + w_edm + w_tlt + w_srn
    w_sar /= w_tot
    w_gps /= w_tot
    w_lev /= w_tot
    w_edm /= w_tot
    w_tlt /= w_tot
    w_srn /= w_tot
    
def read_data():
    global datetime_start
    
    datetime_start = str(datetime.datetime.now())
    print('\nReading data...')
    if(ff_sar[0] !='None'):
        global X_sar, Y_sar, data_sar, err_sar, los_sar,n_sar
        
        n_sar    = []
        X_sar    = []
        Y_sar    = []
        data_sar = []
        err_sar  = []
        los_sar  = []
        n_index  = 0
        
        for i in range(len(ff_sar)):
            ext = ff_sar[i][-3:]
            if(ext == 'txt'):
                try:
                    d_sar = np.loadtxt(ff_sar[i])
                except:
                    d_sar = np.loadtxt(ff_sar[i],skiprows=1)
            if(ext == 'csv'):
                db_sar = pd.read_csv(ff_sar[i])
                d_sar = db_sar.values
            if(ext == 'shp'):
                db_sar = read_shapefile(ff_sar[i])
                d_sar  = db_sar.values
 
            n_sar.append(len(d_sar))
            print('Found ', n_sar[i], 'InSAR data in dataset #',i+1)
            X_sar[n_index:n_sar[i]+n_index]    = d_sar[:,0]
            Y_sar[n_index:n_sar[i]+n_index]    = d_sar[:,1]
            data_sar[n_index:n_sar[i]+n_index] = d_sar[:,2]
            err_sar[n_index:n_sar[i]+n_index]  = d_sar[:,3]
            los_sar[n_index:n_sar[i]+n_index]  = d_sar[:,4:7]
            n_index += n_sar[i]
        X_sar = np.array(X_sar)
        Y_sar = np.array(Y_sar)
        data_sar = np.array(data_sar)
        err_sar = np.array(err_sar)
        los_sar = np.array(los_sar)
        
    if(ff_gps!='None'):
        ext = ff_gps[-3:]
        if(ext == 'txt'):
            try:
                d_gps = np.loadtxt(ff_gps)
            except:
                d_gps = np.loadtxt(ff_gps,skiprows=1)
        if(ext == 'csv'):
            db_gps = pd.read_csv(ff_gps)
            d_gps = db_gps.values
        if(ext == 'shp'):
            db_gps = read_shapefile(ff_gps)
            d_gps  = db_gps.values
            
        global X_gps, Y_gps, data_gps, err_gps
        X_gps    = d_gps[:,0]
        Y_gps    = d_gps[:,1]
        data_gps = d_gps[:,2:5]
        err_gps  = d_gps[:,5:8]
        print('Found ', len(X_gps), 'GPS data')
        
    if(ff_lev!='None'):
        ext = ff_lev[-3:]
        if(ext == 'txt'):
            try:
                d_lev = np.loadtxt(ff_lev)
            except:
                d_lev = np.loadtxt(ff_lev,skiprows=1)
        if(ext == 'csv'):
            db_lev = pd.read_csv(ff_lev)
            d_lev = db_lev.values
        if(ext == 'shp'):
            db_lev = read_shapefile(ff_lev)
            d_lev  = db_lev.values

        global X_lev, Y_lev, data_lev, err_lev
        X_lev    = d_lev[:,0]
        Y_lev    = d_lev[:,1]
        data_lev = d_lev[:,2]
        err_lev  = d_lev[:,3]
        print('Found ', len(X_lev), 'Levelling data')
    
    if(ff_edm!='None'):
        ext = ff_edm[-3:]
        if(ext == 'txt'):
            try:
                d_edm = np.loadtxt(ff_edm)
            except:
                d_edm = np.loadtxt(ff_edm,skiprows=1)
        if(ext == 'csv'):
            db_edm = pd.read_csv(ff_edm)
            d_edm = db_edm.values
        if(ext == 'shp'):
            db_edm = read_shapefile(ff_edm)
            d_edm  = db_edm.values

        global X_edm, Y_edm, data_edm, err_edm
        X_edm      = np.zeros((len(d_edm),2))
        Y_edm      = np.zeros((len(d_edm),2))
        X_edm[:,0] = d_edm[:,0]
        X_edm[:,1] = d_edm[:,2]
        Y_edm[:,0] = d_edm[:,1]
        Y_edm[:,1] = d_edm[:,3]
        data_edm   = d_edm[:,4]
        err_edm    = d_edm[:,5]
        print('Found ', len(X_edm), 'EDM data')

    if(ff_tlt!='None'):
        ext = ff_tlt[-3:]
        if(ext == 'txt'):
            try:
                d_tlt = np.loadtxt(ff_tlt)
            except:
                d_tlt = np.loadtxt(ff_tlt,skiprows=1)
        if(ext == 'csv'):
            db_tlt = pd.read_csv(ff_tlt)
            d_tlt = db_tlt.values
        if(ext == 'shp'):
            db_tlt = read_shapefile(ff_tlt)
            d_tlt  = db_tlt.values

        global X_tlt, Y_tlt, data_tlt, err_tlt
        X_tlt    = d_tlt[:,0]
        Y_tlt    = d_tlt[:,1]
        data_tlt = d_tlt[:,2:4]
        err_tlt  = d_tlt[:,4:6]
        print('Found ', len(X_tlt), 'TILT data')

    if(ff_srn!='None'):
        ext = ff_srn[-3:]
        if(ext == 'txt'):
            try:
                d_srn = np.loadtxt(ff_srn)
            except:
                d_srn = np.loadtxt(ff_srn,skiprows=1)
        if(ext == 'csv'):
            db_srn = pd.read_csv(ff_srn)
            d_srn = db_srn.values
        if(ext == 'shp'):
            db_srn = read_shapefile(ff_srn)
            d_srn  = db_srn.values

        global X_srn, Y_srn, data_srn, err_srn
        X_srn    = d_srn[:,0]
        Y_srn    = d_srn[:,1]
        data_srn = d_srn[:,2]
        err_srn  = d_srn[:,3]
        print('Found ', len(X_srn), 'STRAIN data')
        
def source_info(argument):
    mogi = dict(index = 0,
                name  = 'Mogi Point Source (Mogi, 1958)',
                param_n = 4,
                param_l = ['xcen','ycen','depth','dVol'],
                dpars = None
                )
    mctg = dict(index = 1,
                name  = 'Mogi Finite Volume (McTigue, 1987)',
                param_n = 5,
                param_l = ['xcen','ycen','depth','radius','dP_mu'],
                dpars = None
                )
    fial = dict(index = 2,
                name  = 'Penny-shaped crack (Fialko et al., 2001)',
                param_n = 5,
                param_l = ['xcen','ycen','depth','radius','dP_mu'],
                dpars = None
                )
    yang = dict(index = 3,
                name  = 'Spheroid Finite Volume (Yang et al., 1988)',
                param_n = 8,
                param_l = ['xcen','ycen','depth','s_axis_max','ratio','dP_mu','strike','dip'],
                dpars = None
                )
    davi = dict(index = 4,
                name  = 'Moment Tensor Point Source (Davis, 1986)',
                param_n = 9,
                param_l = ['xcen','ycen','depth','Pxx','Pyy','Pzz','Pxy','Pyz','Pzx'],
                dpars = None
                )
    okad = dict(index = 5,
                name  = 'Rectangular Dislocation (Okada, 1985)',
                param_n = 10,
                param_l = ['xtlc','ytlc','dtlc','length','width','strike','dip','param1','param2','opening'],
                dpars = None
                )
    sources = {
        0: mogi,
        1: mctg,
        2: fial,
        3: yang,
        4: davi,
        5: okad
    }
    source_choice  = sources.get(argument)
    return source_choice
    
def read_VSM_settings(VSM_settings_name):
    global fold_inout
    global ff_gps, ff_lev, ff_edm,ff_tlt,ff_srn
    
    print('\n\n*******************************************************************************')
    print('\n                   VSM exectution begins')
    print('\nStart reading settings of VSM from input file -->')
    print(VSM_settings_name,'\n')
    print('Data -->')
    
    with open(VSM_settings_name, 'r') as reader:
    # Read and print the entire input file line by line
        line = reader.readline()
        fold_inout = line.split()[0]
        
        line = reader.readline()
        try:
            ff_sar[0] = line.split()[0]
            if(ff_sar!='None'): print('InSAR data file #',1,ff_sar[0])
        except:
            pass
        for i in range(1,10):
            try:
                ff_sar.append(line.split()[i])
                print('SAR data file #',i+1,ff_sar[i])
            except:
                break
        
        line = reader.readline()
        try:
            ff_gps = line.split()[0]
            if(ff_gps!='None'): print('GNSS data file', ff_gps)
        except:
            pass
        
        line = reader.readline()
        try:
            ff_lev = line.split()[0]
            if(ff_lev!='None'): print('Levelling data file', ff_lev)
        except:
            pass
        
        line = reader.readline()
        try:
            ff_edm = line.split()[0]
            if(ff_edm!='None'): print('EDM data file', ff_edm)
        except:
            pass
 
        line = reader.readline()
        try:
            ff_tlt = line.split()[0]
            if(ff_tlt!='None'): print('Tilt data file', ff_tlt)
        except:
            pass
        
        line = reader.readline()
        try:
            ff_srn = line.split()[0]
            if(ff_srn!='None'): print('Strain data file', ff_srn)
        except:
            pass
            
        line = reader.readline()
        s = line.split()[0]
        line = reader.readline()
        g = line.split()[0]
        line = reader.readline()
        l = line.split()[0]
        line = reader.readline()
        e = line.split()[0]
        line = reader.readline()
        t = line.split()[0]
        line = reader.readline()
        r = line.split()[0]
        read_weights(s,g,l,e,t,r)
 
        line = reader.readline()
        m = line.split()[0]
        line = reader.readline()
        n = line.split()[0]
        read_elast(m,n)
      
        line = reader.readline()
        ns = int(line.split()[0])

        for i_sources in range(ns):
            line = reader.readline()
            
            option = int(line.split()[0])
            curr_source = source_info(option)
            if(option == 5):
                okada_mode[i_sources] = line.split()[1]
                print('\nSource',i_sources+1,'considered -->\n',curr_source['name'],'with mode',okada_mode[i_sources])
            else:
                print('\nSource',i_sources+1,'considered -->\n',curr_source['name'])
            
            m1=np.zeros(curr_source['param_n'])
            m2=np.zeros(curr_source['param_n'])
            for i_param in range(curr_source['param_n']):
                line = reader.readline()
                m1[i_param] = line.split()[0]
                m2[i_param] = line.split()[1]
            read_params(m1,m2,option)
        
        print('\nThere are', len(parflag),'total free parameters')
        parreal = (np.array(parflag)).sum()
        print('Parameters actually inverted --> ',parreal)
        if(parreal == 1): print('Note that with 1 parameter the 2D marginals won\'t be computed')
        for i in range(len(parflag)):
            if(parflag[i]):
                print(labels[i],bounds[i])
        
        # inversion tool
        line = reader.readline()
        ch   = line.split()[0]
        line = reader.readline()
        sa   = line.split()[0]
        try:
            re   = line.split()[1]
        except:
            re = None
        line = reader.readline()
        it   = line.split()[0]
        read_optim(ch,sa,it,re)
        if(ch_opt==0):
            print('\nInversion tool chosen --> NEIGHBOURHOOD ALGORITHM')
            print('N. of samples: ',samp1, 'N. of re-samples',samp2)
            print('N. of iterations: ',samp3)
            
        if(ch_opt==1):
            print('\nInversion tool chosen --> BAYESIAN INFERENCE')
            print('N. of random walks: ',samp1)
            print('N. of steps: ',samp3)
            if(parreal == 0):
                print('\nUse the Neighbourhhod Algorithm to compute forward models with zero free parameters')
        
        # plot
        try:
            line = reader.readline()
            bo = line.split()[0]
            if(bo.lower() == 'y'): bo = line.split()[1]
        except:
            bo = 2000
        read_plot(bo)
        if(nskip != -1):
            print('\nNumber of burn-in samples for the plots -->',bo,'\n\n')
        else:
            print('\nNo plots this time\n\n')
        
    reader.close()
    
    return

def iVSM():
    read_data()
    print('\nOutput folder is -->\n',fold_inout)
    print('\n                    All input read')
    print('\n*******************************************************************************')
    
    if(ch_opt == 0): my_nei_alg(samp1,samp2,samp3)
    if(ch_opt == 1): my_bay_alg(samp1,samp2,samp3)


# -----------------------------------------------------------------------------
#               DATA INVERSION
# NEIGHBOURHOOD ALGORITHM or BAYESIAN INFERENCE

def synth(params):
    global synth_sar, synth_gps, synth_lev, synth_edm, synth_tlt, synth_srn
        
    if ff_sar[0]!='None': synth_sar = np.zeros( len(X_sar)   )
    if ff_gps   !='None': synth_gps = np.zeros((len(X_gps),3))
    if ff_lev   !='None': synth_lev = np.zeros( len(X_lev)   )
    if ff_edm   !='None': synth_edm = np.zeros( len(X_edm)   )
    if ff_tlt   !='None': synth_tlt = np.zeros((len(X_tlt),2))
    if ff_srn   !='None': synth_srn = np.zeros( len(X_srn)   )
    
    nd = 0
    ndtrue = 0
    for i_sorg in range(len(mysources)):
# set up the parameters for forward functions call
        type_sorg = mysources[i_sorg]['index']
        
        paramsall = np.zeros(mysources[i_sorg]['param_n'])
        
        for iparam in range(mysources[i_sorg]['param_n']):
            if(parflag[nd+iparam]):
                paramsall[iparam] = params[ndtrue]
                ndtrue+=1
            else:
                paramsall[iparam] = bounds[nd+iparam][0]
        
        unknowns = dict(zip(mysources[i_sorg]['param_l'],paramsall))
        #print(unknowns,'\n')
        
        nd += mysources[i_sorg]['param_n']
        unknowns['nu'] = ni

        if(type_sorg == 5):
            unknowns['opt'] = okada_mode[i_sorg]

# separate calls for each dataset type        
        
        if ff_sar[0]!='None':
            if(type_sorg == 0):
                ux,uy,uz = forward.mogi(   X_sar, Y_sar, **unknowns)
            if(type_sorg == 1):
                ux,uy,uz = forward.mctigue(X_sar, Y_sar, **unknowns)
            if(type_sorg == 2):
                ux,uy,uz = forward.fialko( X_sar, Y_sar, **unknowns)
            if(type_sorg == 3):
                ux,uy,uz = forward.yang(   X_sar, Y_sar, **unknowns)
            if(type_sorg == 4):
                ux,uy,uz = forward.davis(  X_sar, Y_sar, **unknowns)
            if(type_sorg == 5):
                ux,uy,uz = forward.okada(  X_sar, Y_sar, **unknowns)
               
            utot = np.array([ux,uy,uz])
            utot = np.transpose(utot)
            ulos = np.sum(utot*los_sar,axis=1)
            synth_sar += ulos

        if ff_gps!='None':
            if(type_sorg == 0):
                ux,uy,uz = forward.mogi(   X_gps, Y_gps, **unknowns)
            if(type_sorg == 1):
                ux,uy,uz = forward.mctigue(X_gps, Y_gps, **unknowns)
            if(type_sorg == 2):
                ux,uy,uz = forward.fialko( X_gps, Y_gps, **unknowns)
            if(type_sorg == 3):
                ux,uy,uz = forward.yang(   X_gps, Y_gps, **unknowns)
            if(type_sorg == 4):
                ux,uy,uz = forward.davis(  X_gps, Y_gps, **unknowns)
            if(type_sorg == 5):
                ux,uy,uz = forward.okada(  X_gps, Y_gps, **unknowns)
                
            synth_gps[:,0] += ux
            synth_gps[:,1] += uy
            synth_gps[:,2] += uz
                      
        if ff_lev!='None':
            if(type_sorg == 0):
                ux,uy,uz = forward.mogi(   X_lev, Y_lev, **unknowns)
            if(type_sorg == 1):
                ux,uy,uz = forward.mctigue(X_lev, Y_lev, **unknowns)
            if(type_sorg == 2):
                ux,uy,uz = forward.fialko( X_lev, Y_lev, **unknowns)
            if(type_sorg == 3):
                ux,uy,uz = forward.yang(   X_lev, Y_lev, **unknowns)
            if(type_sorg == 4):
                ux,uy,uz = forward.davis(  X_lev, Y_lev, **unknowns)
            if(type_sorg == 5):
                ux,uy,uz = forward.okada(  X_lev, Y_lev, **unknowns)
            
            synth_lev += uz

        if ff_edm!='None':
            X_edm0 = X_edm[:,0]
            Y_edm0 = Y_edm[:,0]
            X_edm1 = X_edm[:,1]
            Y_edm1 = Y_edm[:,1]
            if(type_sorg == 0):
                ux0,uy0,uz0 = forward.mogi(   X_edm0, Y_edm0, **unknowns)
                ux1,uy1,uz1 = forward.mogi(   X_edm1, Y_edm1, **unknowns)
            if(type_sorg == 1):
                ux0,uy0,uz0 = forward.mctigue(X_edm0, Y_edm0, **unknowns)
                ux1,uy1,uz1 = forward.mctigue(X_edm1, Y_edm1, **unknowns)
            if(type_sorg == 2):
                ux0,uy0,uz0 = forward.fialko( X_edm0, Y_edm0, **unknowns)
                ux1,uy1,uz1 = forward.fialko( X_edm1, Y_edm1, **unknowns)
            if(type_sorg == 3):
                ux0,uy0,uz0 = forward.yang(   X_edm0, Y_edm0, **unknowns)
                ux1,uy1,uz1 = forward.yang(   X_edm1, Y_edm1, **unknowns)
            if(type_sorg == 4):
                ux0,uy0,uz0 = forward.davis(  X_edm0, Y_edm0, **unknowns)
                ux1,uy1,uz1 = forward.davis(  X_edm1, Y_edm1, **unknowns)
            if(type_sorg == 5):
                ux0,uy0,uz0 = forward.okada(  X_edm0, Y_edm0, **unknowns)
                ux1,uy1,uz1 = forward.okada(  X_edm1, Y_edm1, **unknowns)
                
            ddxx  = X_edm0 - X_edm1
            ddyy  = Y_edm0 - Y_edm1
            dist1 = np.sqrt(ddxx**2 + ddyy**2)
            ddxx  += ux0 - ux1
            ddyy  += uy0 - uy1
            dist2 = np.sqrt(ddxx**2 + ddyy**2)
            synth_edm += (dist2 - dist1)

        if ff_tlt!='None':
            delta = 2.
            if(type_sorg == 0):
                uxpX,uypX,uzpX = forward.mogi(   X_tlt+delta, Y_tlt, **unknowns)
                uxmX,uymX,uzmX = forward.mogi(   X_tlt-delta, Y_tlt, **unknowns)
                uxpY,uypY,uzpY = forward.mogi(   X_tlt, Y_tlt+delta, **unknowns)
                uxmY,uymY,uzmY = forward.mogi(   X_tlt, Y_tlt-delta, **unknowns)
            if(type_sorg == 1):
                uxpX,uypX,uzpX = forward.mctigue(X_tlt+delta, Y_tlt, **unknowns)
                uxmX,uymX,uzmX = forward.mctigue(X_tlt-delta, Y_tlt, **unknowns)
                uxpY,uypY,uzpY = forward.mctigue(X_tlt, Y_tlt+delta, **unknowns)
                uxmY,uymY,uzmY = forward.mctigue(X_tlt, Y_tlt-delta, **unknowns)
            if(type_sorg == 2):
                uxpX,uypX,uzpX = forward.fialko( X_tlt+delta, Y_tlt, **unknowns)
                uxmX,uymX,uzmX = forward.fialko( X_tlt-delta, Y_tlt, **unknowns)
                uxpY,uypY,uzpY = forward.fialko( X_tlt, Y_tlt+delta, **unknowns)
                uxmY,uymY,uzmY = forward.fialko( X_tlt, Y_tlt-delta, **unknowns)
            if(type_sorg == 3):
                uxpX,uypX,uzpX = forward.yang(   X_tlt+delta, Y_tlt, **unknowns)
                uxmX,uymX,uzmX = forward.yang(   X_tlt-delta, Y_tlt, **unknowns)
                uxpY,uypY,uzpY = forward.yang(   X_tlt, Y_tlt+delta, **unknowns)
                uxmY,uymY,uzmY = forward.yang(   X_tlt, Y_tlt-delta, **unknowns)
            if(type_sorg == 4):
                uxpX,uypX,uzpX = forward.davis(  X_tlt+delta, Y_tlt, **unknowns)
                uxmX,uymX,uzmX = forward.davis(  X_tlt-delta, Y_tlt, **unknowns)
                uxpY,uypY,uzpY = forward.davis(  X_tlt, Y_tlt+delta, **unknowns)
                uxmY,uymY,uzmY = forward.davis(  X_tlt, Y_tlt-delta, **unknowns)
            if(type_sorg == 5):
                uxpX,uypX,uzpX = forward.okada(  X_tlt+delta, Y_tlt, **unknowns)
                uxmX,uymX,uzmX = forward.okada(  X_tlt-delta, Y_tlt, **unknowns)
                uxpY,uypY,uzpY = forward.okada(  X_tlt, Y_tlt+delta, **unknowns)
                uxmY,uymY,uzmY = forward.okada(  X_tlt, Y_tlt-delta, **unknowns)

            synth_tlt[:,0] += (uzpX-uzmX)*1.e6/delta/2.
            synth_tlt[:,1] += (uzpY-uzmY)*1.e6/delta/2.

        if ff_srn!='None':
            delta = 2.
            fact = (1.-2.*ni)/(1.+ni)
            if(type_sorg == 0):
                uxpX,uypX,uzpX = forward.mogi(   X_srn+delta, Y_srn, **unknowns)
                uxmX,uymX,uzmX = forward.mogi(   X_srn-delta, Y_srn, **unknowns)
                uxpY,uypY,uzpY = forward.mogi(   X_srn, Y_srn+delta, **unknowns)
                uxmY,uymY,uzmY = forward.mogi(   X_srn, Y_srn-delta, **unknowns)
            if(type_sorg == 1):
                uxpX,uypX,uzpX = forward.mctigue(X_srn+delta, Y_srn, **unknowns)
                uxmX,uymX,uzmX = forward.mctigue(X_srn-delta, Y_srn, **unknowns)
                uxpY,uypY,uzpY = forward.mctigue(X_srn, Y_srn+delta, **unknowns)
                uxmY,uymY,uzmY = forward.mctigue(X_srn, Y_srn-delta, **unknowns)
            if(type_sorg == 2):
                uxpX,uypX,uzpX = forward.fialko( X_srn+delta, Y_srn, **unknowns)
                uxmX,uymX,uzmX = forward.fialko( X_srn-delta, Y_srn, **unknowns)
                uxpY,uypY,uzpY = forward.fialko( X_srn, Y_srn+delta, **unknowns)
                uxmY,uymY,uzmY = forward.fialko( X_srn, Y_srn-delta, **unknowns)
            if(type_sorg == 3):
                uxpX,uypX,uzpX = forward.yang(   X_srn+delta, Y_srn, **unknowns)
                uxmX,uymX,uzmX = forward.yang(   X_srn-delta, Y_srn, **unknowns)
                uxpY,uypY,uzpY = forward.yang(   X_srn, Y_srn+delta, **unknowns)
                uxmY,uymY,uzmY = forward.yang(   X_srn, Y_srn-delta, **unknowns)
            if(type_sorg == 4):
                uxpX,uypX,uzpX = forward.davis(  X_srn+delta, Y_srn, **unknowns)
                uxmX,uymX,uzmX = forward.davis(  X_srn-delta, Y_srn, **unknowns)
                uxpY,uypY,uzpY = forward.davis(  X_srn, Y_srn+delta, **unknowns)
                uxmY,uymY,uzmY = forward.davis(  X_srn, Y_srn-delta, **unknowns)
            if(type_sorg == 5):
                uxpX,uypX,uzpX = forward.okada(  X_srn+delta, Y_srn, **unknowns)
                uxmX,uymX,uzmX = forward.okada(  X_srn-delta, Y_srn, **unknowns)
                uxpY,uypY,uzpY = forward.okada(  X_srn, Y_srn+delta, **unknowns)
                uxmY,uymY,uzmY = forward.okada(  X_srn, Y_srn-delta, **unknowns)

            synth_srn += fact*(uxpX - uxmX + uypY - uymY)*1.e6/delta/2.
                
    return

# -----------------------------------------------------------------------------
# NEIGHBOURHOOD ALGORITHM

def chi_square(dd,tt,ee):
    chi2 = np.sum(((dd-tt)/ee)**2)
    chi2 /= dd.size
    return chi2

def na_obj(params):
    func = None
    
    wchi2_sar = 0.
    wchi2_gps = 0.
    wchi2_lev = 0.
    wchi2_edm = 0.
    wchi2_tlt = 0.
    wchi2_srn = 0.
    
    synth(params)
             
    if ff_sar[0]!='None':
        chi2_sar = chi_square(data_sar,synth_sar,err_sar)
        wchi2_sar = w_sar*chi2_sar
    if ff_gps!='None':   
        chi2_gps = chi_square(data_gps,synth_gps,err_gps)
        wchi2_gps = w_gps*chi2_gps
    if ff_lev!='None':   
        chi2_lev = chi_square(data_lev,synth_lev,err_lev)
        wchi2_lev = w_lev*chi2_lev
    if ff_edm!='None':   
        chi2_edm = chi_square(data_edm,synth_edm,err_edm)
        wchi2_edm = w_edm*chi2_edm
    if ff_tlt!='None':   
        chi2_tlt = chi_square(data_tlt,synth_tlt,err_tlt)
        wchi2_tlt = w_tlt*chi2_tlt
    if ff_srn!='None':   
        chi2_srn = chi_square(data_srn,synth_srn,err_srn)
        wchi2_srn = w_srn*chi2_srn
        
    func = wchi2_sar + wchi2_gps + wchi2_lev + wchi2_edm + wchi2_tlt + wchi2_srn
    
    return func

def my_nei_alg(samp1,samp2,samp3):
    print('\nNEIGHBOURHOOD ALGORITM running...')
    num_samp   = samp1
    num_resamp = samp2
    num_iter   = samp3
    
    ndim = len(bounds_true)
    
    srch_NA = NAsearch.Searcher(
        objective   = na_obj,
        limits      = bounds_true,
        num_samp    = num_samp,
        num_resamp  = num_resamp,
        names       = labels_true,
        maximize    = False,
        verbose     = True
        )

    srch_NA.update(num_iter)
    flat_samples = srch_NA.sample_dataframe
    flat_samples = flat_samples.to_numpy()
    
    best = flat_samples[0,2:]
    print('\nBEST PARAMETERS\n',best)
    
### Write results of NA    
    write_results(flat_samples[:,2:],
                  truths = best,
                  extra   = flat_samples[:,:2])
    
    if(plotyn == True):

### Plot parameters sampling
        if(ndim >= 1):
            plot_params(flat_samples,
                    labels   = labels_true,
                    filename = os.path.join(fold_inout,'VSM_params.png'))

### Plot 1D 2D Density Distributions
        if(ndim >= 2):
            plot_1D_2D(flat_samples[:,2:],
                   truths   = best,
                   labels   = labels_true,
                   filename = os.path.join(fold_inout,'VSM_1D2D.png'), 
                   skip     = nskip)

    return flat_samples

# -----------------------------------------------------------------------------
# BAYESIAN INFERENCE

def log_likelihood(dd,tt,ee):
    chi2 = chi_square(dd,tt,ee)
    ee2 = ee**2
    log_lik = -0.5 * (dd.size* chi2 + np.sum(np.log(ee2*2.*np.pi)) )
    
    return log_lik
 
def log_likelitot(params):
    
    log_lik_sar = 0.
    log_lik_gps = 0.
    log_lik_lev = 0.
    log_lik_edm = 0.
    log_lik_tlt = 0.
    log_lik_srn = 0.
    
    synth(params)
    
    if ff_sar[0]!='None':   
        log_lik_sar = log_likelihood(data_sar,synth_sar,err_sar)
    if ff_gps!='None': 
        log_lik_gps = log_likelihood(data_gps,synth_gps,err_gps)
    if ff_lev!='None':   
        log_lik_lev = log_likelihood(data_lev,synth_lev,err_lev)
    if ff_edm!='None':   
        log_lik_edm = log_likelihood(data_edm,synth_edm,err_edm)
    if ff_tlt!='None':   
        log_lik_tlt = log_likelihood(data_tlt,synth_tlt,err_tlt)
    if ff_srn!='None':   
        log_lik_srn = log_likelihood(data_srn,synth_srn,err_srn)
        
    log_lik_tot = log_lik_sar*w_sar + log_lik_gps*w_gps + log_lik_lev*w_lev + \
                  log_lik_edm*w_edm + log_lik_tlt*w_tlt + log_lik_srn*w_srn
                  
    return log_lik_tot

def log_prior(params):
    k = 0
    for i in range(len(bounds_true)):
        if bounds_true[i][0] <= params[i] <= bounds_true[i][1]:
            k+=1
    if(k != len(bounds_true)):
        ret = -np.inf
    else:
        ret = 0
    
    return ret

def log_probability(params):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelitot(params)

def my_bay_alg(samp1,samp2,samp3,excl=80,nbins=20):
    
    nwalkers = samp1
    thin     = samp2
    nsteps   = samp3
    
    print('\nBAYESIAN INFERENCE running...')
    ndim = len(bounds_true)
    
    pos0 = np.random.random((nwalkers, ndim))
    pos = np.zeros((nwalkers,ndim))
    for i in range(ndim):
        pos[:,i] = pos0[:,i]*(bounds_true[i][1]-bounds_true[i][0])+bounds_true[i][0]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
  
    sampler.run_mcmc(pos,nsteps,progress=True)

    flat_samples = sampler.get_chain(discard=excl, thin=thin, flat=True)

    mean  = np.zeros(ndim)
    sigma = np.zeros(ndim)
    print('\nMEAN PARAMETERS and associated STANDARD DEVIATION')
    for i in range(ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        mean[i] = mcmc[1]
        q = np.diff(mcmc)
        sigma[i] = np.sum(q)/2.
        print(labels_true[i],mean[i],'+-',sigma[i])

### Write results of BO 
    write_results(flat_samples,
                  truths = mean,
                  sigma   = sigma)

    if(plotyn == True):
        
### Plot parameters sampling
        if(ndim >= 1):
            plot_params(flat_samples,
                    labels   = labels_true,
                    filename = os.path.join(fold_inout,'VSM_params.png'))
        
### Plot 1D 2D Posterior Probability Density Distribution    
        if(ndim >= 2):
            plot_1D_2D(flat_samples,
                   truths   = mean,
                   labels   = labels_true,
                   filename = os.path.join(fold_inout,'VSM_1D2D.png'), 
                   skip     = nskip)

    return flat_samples

# -----------------------------------------------------------------------------
# WRITE RESULTS
def write_results(samples, truths, sigma=None, extra=None):
    
    ndim = len(truths)

    nbins = 20
    
### Reconstruct best fit model with synthetic data associated
    chi2_final = na_obj(truths)
    print('\nTotal misfit of the optimal model {:10.4f}'.format(chi2_final))

### Save best fit (NA) or mean model with 1-sigma (BI)
    if(ch_opt == 0 and ndim != 0):
        filename_temp = os.path.join(fold_inout,'VSM_best.csv')
        hd=''
        for i in range(ndim):
            hd = hd+labels_true[i]+','
            best = np.array([truths])
        hd = hd[:-1]
        np.savetxt(filename_temp,best,delimiter=',',header=hd,comments='',fmt='%.7e')  
    if(ch_opt == 1 and ndim != 0):
        filename_temp = os.path.join(fold_inout,'VSM_mean.csv')
        hd = ''
        all_mean=[]
        for i in range(ndim):
            all_mean += truths[i],sigma[i]
            hd += str(labels_true[i])+','+str(labels_true[i])+'_sigma,'
        hd = hd[:-1]
        all_mean = np.array([all_mean])
        np.savetxt(filename_temp,all_mean,header=hd,comments='', fmt = (((2*ndim)-1)*('%.7e,')+'%.7e'))

### Save all models generated
    if(ch_opt == 0 and ndim != 0):
        filename_temp = os.path.join(fold_inout,'VSM_models.csv')
        dff = np.concatenate((extra,samples),axis=1)
        hd='func,iter,'+hd
        np.savetxt(filename_temp,dff,header=hd,comments='',fmt='%.6e'+',%i'+',%.7e'*ndim)
    if(ch_opt == 1 and ndim != 0):
        hd=''
        for i in range(ndim):
            hd = hd+labels_true[i]+','
            best = np.array([truths])
        hd = hd[:-1]
        filename_temp = os.path.join(fold_inout,'VSM_models.csv')
        np.savetxt(filename_temp,samples,header=hd,delimiter=',',comments='',fmt='%.7e')

### Save all models generated 1D marginal PPD 
    filename_temp = os.path.join(fold_inout,'VSM_1D.csv')
    h1d   = np.zeros((nbins,ndim))
    e1d   = np.zeros((nbins+1,ndim))
    for i in range(ndim):
        h1d[:,i], e1d[:,i] = np.histogram(samples[:,i], bins= nbins, density=True)
    
    hd = ''
    all_1D = np.zeros((nbins,ndim*2))
    for i in range(ndim):
        edge_bin = e1d[1,i] - e1d[0,i]
        edges = e1d[:-1,i] + edge_bin
        all_1D[:,i*2]= edges
        all_1D[:,i*2+1] = h1d[:,i]*edge_bin # multiply for edge bin to normalize to parameter range
        hd += str(labels_true[i])+','+str(labels_true[i])+'_1D_PPD,'
    hd = hd[:-1]
    np.savetxt(filename_temp,all_1D,header=hd,comments='', fmt = (((2*ndim)-1)*('%.5e,')+'%.5e'))

### Save all models generated 2D marginal PPD
    if(ndim >= 2):
        filename_temp = os.path.join(fold_inout,'VSM_2D.txt')
       
        with open(filename_temp, "w") as f:
            for i in range(ndim):
                for j in range(i+1,ndim):
                    h2d, e2dx,e2dy = np.histogram2d(samples[:,i],samples[:,j],bins=nbins,density=True)
                    xm = e2dx.min()
                    xp = e2dx.max()
                    ym = e2dy.min()
                    yp = e2dy.max()
                
                    hd= '2D Marginal for parameters {} and {} with ranges {:.4e} {:.4e} and {:.4e} {:.4e}'.format(labels_true[i],labels_true[j],xm,xp,ym,yp)
                    for l in range(nbins):
                        line = h2d[l]
                        line = np.array([line])
                        if(l == 0):
                            np.savetxt(f, line,header=hd, fmt='%.3e')
                        else:
                            np.savetxt(f, line, fmt='%.3e')

### Save synthetic data
    if ff_sar[0] != 'None':
        n_index = 0
        for i in range(len(n_sar)):
            if len(n_sar) == 1:
                filename_temp = os.path.join(fold_inout,'VSM_synth_sar.csv')
                all_sar = np.asarray([X_sar,Y_sar,synth_sar,data_sar,err_sar,los_sar[:,0],los_sar[:,1],los_sar[:,2]])
            else:
                filename_temp = os.path.join(fold_inout,'VSM_synth_sar'+str(i+1)+'.csv')
                xx = X_sar[n_index:n_sar[i]+n_index]
                yy = Y_sar[n_index:n_sar[i]+n_index]
                dd = data_sar[n_index:n_sar[i]+n_index]
                ee = err_sar[n_index:n_sar[i]+n_index]
                ss = synth_sar[n_index:n_sar[i]+n_index]
                le = los_sar[n_index:n_sar[i]+n_index,0]
                ln = los_sar[n_index:n_sar[i]+n_index,1]
                lz = los_sar[n_index:n_sar[i]+n_index,2]
                n_index += n_sar[i]
                all_sar = np.asarray([xx,yy,ss,dd,ee,le,ln,lz])
            hd='coo_X,coo_Y,synth_sar,data_sar,err_sar,los_E,los_N,los_Z'
            form  = 2*'%.2f,'+3*'%.5e,'+3*'%.3f,'
            form= form[:-1]
            np.savetxt(filename_temp,all_sar.transpose(),delimiter=',',header=hd,comments='',fmt=form)   
    if ff_gps != 'None':
        filename_temp = os.path.join(fold_inout,'VSM_synth_gps.csv')
        all_gps = np.asarray([X_gps,Y_gps,synth_gps[:,0],synth_gps[:,1],synth_gps[:,2],
                                           data_gps[:,0], data_gps[:,1], data_gps[:,2],
                                            err_gps[:,0],  err_gps[:,1],  err_gps[:,2]])
        hd='coo_X,coo_Y,synth_gps_E,synth_gps_N,synth_gps_Z,data_gps_E,data_gps_N,data_gps_Z,err_gps_E,err_gps_N,err_gps_Z'
        form  = 2*'%.2f,'+9*'%.5e,'
        form= form[:-1]
        np.savetxt(filename_temp,all_gps.transpose(),delimiter=',',header=hd,comments='',fmt=form)
    if ff_lev != 'None':
        filename_temp = os.path.join(fold_inout,'VSM_synth_lev.csv')
        all_lev = np.asarray([X_lev,Y_lev,synth_lev,data_lev,err_lev])
        hd='coo_X,coo_Y,synth_lev,data_lev,err_lev'
        form  = 2*'%.2f,'+3*'%.5e,'
        form= form[:-1]
        np.savetxt(filename_temp,all_lev.transpose(),delimiter=',',header=hd,comments='',fmt=form)
    if ff_edm != 'None':
        filename_temp = os.path.join(fold_inout,'VSM_synth_edm.csv')
        all_edm = np.asarray([X_edm[:,0],Y_edm[:,0],X_edm[:,1],Y_edm[:,1],synth_edm,data_edm,err_edm])
        hd='coo_X0,coo_Y0,coo_X1,coo_Y1,synth_sar,data_sar,err_sar'
        form  = 4*'%.2f,'+3*'%.5e,'
        form= form[:-1]
        np.savetxt(filename_temp,all_edm.transpose(),delimiter=',',header=hd,comments='',fmt=form)
    if ff_tlt != 'None':
        filename_temp = os.path.join(fold_inout,'VSM_synth_tlt.csv')
        all_tlt = np.asarray([X_tlt,Y_tlt,synth_tlt[:,0],synth_tlt[:,1],
                                           data_tlt[:,0], data_tlt[:,1],
                                            err_tlt[:,0],  err_tlt[:,1]])
        hd='coo_X,coo_Y,synth_tlt_1,synth_tlt_2,data_tlt_1,data_tlt_2,err_tlt_1,err_tlt_2'
        form  = 2*'%.2f,'+6*'%.5e,'
        form= form[:-1]
        np.savetxt(filename_temp,all_tlt.transpose(),delimiter=',',header=hd,comments='',fmt=form)   
    if ff_srn != 'None':
        filename_temp = os.path.join(fold_inout,'VSM_synth_srn.csv')
        all_srn = np.asarray([X_srn,Y_srn,synth_srn,data_srn,err_srn])
        hd='coo_X,coo_Y,synth_srn,data_srn,err_srn'
        form  = 2*'%.2f,'+3*'%.5e,'
        form= form[:-1]
        np.savetxt(filename_temp,all_srn.transpose(),delimiter=',',header=hd,comments='',fmt=form)   
        
### Save all info from the current run of VSM
    filename_temp = os.path.join(fold_inout,'VSM.log')
    stars = '*******************************************************************************\n'
    f = open(filename_temp, "w")
    #f.writelines([l for l in open("./license.lic").readlines(2100)])
    f.write(stars)
    f.write('   VSM\n\n')
    f.write('   Volcanic and Seismic source Modelling\n\n')
    f.write('   Author: Elisa Trasatti, elisa.trasatti@ingv.it\n')
    f.write('   Istituto Nazionale di Geofisica e Vulcanologia - Rome (Italy)\n\n')
    f.write('   License: E. Trasatti, covered by GNU-GPL License\n\n')
    f.write('   Redistribution and use in source and binary forms, with or without\n')
    f.write('   modification, are permitted provided that the following conditions are met:\n')
    f.write('      * Redistributions of source code must retain the above copyright \n')
    f.write('        notice, this list of conditions and the following disclaimer.\n')
    f.write('      * Redistributions in binary form must reproduce the above copyright\n')
    f.write('        notice, this list of conditions and the following disclaimer in\n')
    f.write('        the documentation and/or other materials provided with the distribution.\n\n')
    f.write('   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" \n')
    f.write('   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n') 
    f.write('   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE \n')
    f.write('   ARE DISCLAIMED. IN NO EVENT SHALL THE COPY RIGHT OWNER OR CONTRIBUTORS BE\n')
    f.write('   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR \n')
    f.write('   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF \n')
    f.write('   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS\n')
    f.write('   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN \n')
    f.write('   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) \n')
    f.write('   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE\n')
    f.write('   POSSIBILITY OF SUCH DAMAGE.\n')
    f.write(stars)
    ##
    f.write('\n\n'+stars)
    f.write('VSM launched on '+datetime_start+'\n')
    f.write('VSM ended on    '+str(datetime.datetime.now())+'\n')
    f.write(stars+'\n')
    ##
    f.write('\n'+stars)
    f.write('  Misfits \n')
    f.write(stars)
    f.write('Dataset\tn.Data\tWeight\tMisfit_null\tMisfit\n')
    if ff_sar[0] != 'None':
        ll = len(X_sar)
        chi2_null = chi_square(data_sar,0,err_sar)
        chi2_sar  = chi_square(data_sar,synth_sar,err_sar)
        f.write('SAR\t'+str(ll)+'\t%8.4f\t%8.3f\t%8.3f \n' % (w_sar,chi2_null, chi2_sar))
    else:
        f.write('SAR\t - \t -\t    -\t\t  -\n')
    if ff_gps != 'None':
        ll = len(X_gps)
        chi2_null = chi_square(data_gps,0,err_gps)
        chi2_gps  = chi_square(data_gps,synth_gps,err_gps)
        f.write('GPS\t'+str(ll)+'\t%8.4f\t%8.3f\t%8.3f \n' % (w_gps,chi2_null, chi2_gps))
    else:
        f.write('GPS\t - \t -\t    -\t\t  -\n')
    if ff_lev != 'None':
        ll = len(X_lev)
        chi2_null = chi_square(data_lev,0,err_lev)
        chi2_lev  = chi_square(data_lev,synth_lev,err_lev)
        f.write('LEV\t'+str(ll)+'\t%8.4f\t%8.3f\t%8.3f \n' % (w_lev,chi2_null, chi2_lev))
    else:
        f.write('LEV\t - \t -\t    -\t\t  -\n')
    if ff_edm != 'None':
        ll = len(X_edm)
        chi2_null = chi_square(data_edm,0,err_edm)
        chi2_edm  = chi_square(data_edm,synth_edm,err_edm)
        f.write('EDM\t'+str(ll)+'\t%8.4f\t%8.3f\t%8.3f \n' % (w_edm,chi2_null, chi2_edm))
    else:
        f.write('EDM\t - \t -\t    -\t\t  -\n')
    if ff_tlt != 'None':
        ll = len(X_tlt)
        chi2_null = chi_square(data_tlt,0,err_tlt)
        chi2_tlt  = chi_square(data_tlt,synth_tlt,err_tlt)
        f.write('TLT\t'+str(ll)+'\t%8.4f\t%8.3f\t%8.3f \n' % (w_tlt,chi2_null, chi2_tlt))
    else:
        f.write('TLT\t - \t -\t    -\t\t  -\n')
    if ff_srn != 'None':
        ll = len(X_srn)
        chi2_null = chi_square(data_srn,0,err_srn)
        chi2_srn  = chi_square(data_srn,synth_srn,err_srn)
        f.write('SRN\t'+str(ll)+'\t%8.4f\t%8.3f\t%8.3f \n' % (w_srn,chi2_null, chi2_srn))
    else:
        f.write('SRN\t - \t -\t    -\t\t  -\n')
        
    f.write('\nTotal misfit of the preferred model '+'%8.3f' % (chi2_final))
    ##    
    f.write('\n\n'+stars)
    f.write('  Forward Models \n')
    f.write(stars)
    for i_sorg in range(len(mysources)):
        if(mysources[i_sorg]['index']!= 5):
            f.write(str(mysources[i_sorg]['name'])+'\n')
        else:
            f.write(str(mysources[i_sorg]['name'])+' with mode '+okada_mode[i_sorg]+'\n')
    ##    
    f.write('\n'+stars)
    f.write('  Unknowns \n')
    f.write(stars)
    if(ch_opt==0): f.write('Parameter\tBest\t\tRange\t\t\tStatus\n')
    if(ch_opt==1): f.write('Parameter\tMean\t\tSigma\t\tRange\t\t\tStatus\n')
       
    ndtrue = 0
    for i in range(len(labels)):
        if(len(labels[i]) < 8):
            space='\t\t'
        else:
            space='\t'
        if(parflag[i]):
            if(ch_opt == 0):
                row = labels[i]+space+'%.5e\t' % (truths[ndtrue])+str(bounds[i])+'\t true\n'
            if(ch_opt == 1):
                row = labels[i]+space+'%.5e\t%.5e\t' % (truths[ndtrue],sigma[ndtrue])+str(bounds[i])+'\t true\n'
            ndtrue+=1
        else:
            if(ch_opt == 0):
                row = labels[i]+space+'%.5e\t' % (bounds[i][0])+str(bounds[i])+' false\n'
            if(ch_opt == 1):
                row = labels[i]+space+'%.5e\t0.0\t\t' % (bounds[i][0])+str(bounds[i])+' false\n'
        f.write(row)
        
    ##    
    f.write('\n'+stars)
    f.write('  Data Inversion \n')
    f.write(stars)
    if(ch_opt==0):
        f.write('Neighbourhood Algorithm \n')
        f.write('n.Samples\tn.Resamples\tn.Iterations\n')
        f.write(str(samp1)+'\t\t'+str(samp2)+'\t\t'+str(samp3))
        f.write('\n\nTotal number of models generated '+str(samp1*samp3))
    if(ch_opt==1):
        f.write('Bayesian Inference \n')
        f.write('n.Rand.Walks\tn.Steps\n')
        f.write(str(samp1)+'\t\t'+str(samp3))
        f.write('\n\nTotal number of models generated '+str(samp1*samp3))
    ##    
    f.write('\n\n'+stars)
    f.write('  All input - VSM style \n')
    f.write(stars)
    f.write(fold_inout+'\n')
    
    try:
        if len(ff_sar)<2:
            f.write(ff_sar[0]+'\n')
        else:
            ff_sar_temp = ''
            for i in range(len(ff_sar)):
                ff_sar_temp += ff_sar[i]+' '
            f.write(ff_sar_temp+'\n')
    except:
        f.write(ff_sar[0]+'\n')
    
    f.write(ff_gps+'\n'+ff_lev+'\n'+ff_edm+'\n'+ff_tlt+'\n'+ff_srn+'\n')
    f.write(str(w_sar)+'\n'+str(w_gps)+'\n'+str(w_lev)+'\n'+str(w_edm)+'\n'+str(w_tlt)+'\n'+str(w_srn)+'\n')
    f.write(('%.2e\n' % mu)+str(ni)+'\n')
    f.write(str(len(mysources))+'\n')
    kk = 0
    for i_sorg in range(len(mysources)):
        if(mysources[i_sorg]['index']!= 5):
            f.write(str(mysources[i_sorg]['index'])+'\n')
        else:
            f.write(str(mysources[i_sorg]['index'])+' '+okada_mode[i_sorg]+'\n')
       
        for k in range(mysources[i_sorg]['param_n']):
            row = str(list(bounds[kk])[0])+'\t'+str(list(bounds[kk])[1])+'\t'+mysources[i_sorg]['param_l'][k]+'\n'
            kk +=1            
            f.write(row)
    f.write(str(ch_opt)+'\n')
    f.write(str(samp1)+' '+str(samp2)+'\n')
    f.write(str(samp3)+'\n')
    f.write(str(nskip)+'\n')
    
    f.close()
     
# -----------------------------------------------------------------------------
# PLOTS 

def plot_1D_2D(what, truths, labels,filename = None, skip = 50):

        print('\nPlotting 1D and 2D parameter marginals \n')
        
        fig2 = plt.figure(figsize=(14,14))
        corner.corner(what[:-skip],truths=list(truths),fig= fig2,labels = labels)
        
        if filename:
            plt.savefig(filename)
        else:
            fig2.show()

def plot_params(what,labels,filename=None):
    
        print('\nPlotting parameters sampling \n')

        ndim = len(labels)
        fig3, axes = plt.subplots(ndim, figsize=(14, 14), sharex=True)
        
        if(ch_opt==0):
            samples_unsort = what[np.argsort(what[:, 1])]
            for i in range(ndim):
                if(ndim >= 2):
                    ax = axes[i]
                else:
                    ax = axes
                ax.plot(samples_unsort[:, i+2], "k", alpha=0.3)
                ax.set_xlim(0, len(samples_unsort)-2)
                ax.set_ylabel(labels[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)
            if(ndim >= 2):
                axes[-1].set_xlabel("Samples")
            else:
                axes.set_xlabel("Samples")
            
        if(ch_opt==1):
            for i in range(ndim):
                if(ndim >= 2):
                    ax = axes[i]
                else:
                    ax = axes
                ax.plot(what[:, i], "k", alpha=0.3)
                ax.set_xlim(0, len(what))
                ax.set_ylabel(labels[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)
            if(ndim >= 2):
                axes[-1].set_xlabel("Samples")
            else:
                axes.set_xlabel("Samples")

        if filename:
            plt.savefig(filename)
        else:
            fig3.show()
                
# -----------------------------------------------------------------------------
# MAIN            
if __name__ == "__main__":

    filename_in = '../USECASE/mogi1NA/VSM_input_mogi.txt'
    
    read_VSM_settings(filename_in)
    iVSM()
