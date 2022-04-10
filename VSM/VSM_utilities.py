#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VSM - Version 1.0

   Volcanic and Seismic source Modelling

   Author:  Elisa Trasatti, elisa.trasatti@ingv.it
   Istituto Nazionale di Geofisica e Vulcanologia - Rome (Italy)

   Last update:  April 2022

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
   
*******************************************************************************

   BIBLIOGRAPHY

   Amoruso, A., Crescentini, L. (2013). Analytical models of volcanic 
       ellipsoidal expansion sources. Ann. Geophys.
       https://doi.org/10.4401/ag-6441
   Battaglia, M., Cervelli, P.F., Murray, J.R. (2013). Modeling Crustal 
       Deformation near Active Faults and Volcanic Centers—A Catalog of 
       Deformation Models, Book 13, Volcanic Monitoring, USGS.
   Okada, Y. (1985). Surface deformation due to shear and tensile faults in a 
       half-space. B. Seismol. Soc. Am. 75, 1135-1154. 
       
*******************************************************************************

    This file contains several utilities to perform post-processing.
    Uncomment to use one of them, or use VSM_utilities.ipynb.
    
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

#import sys
#sys.path.append('path-to-VSM-folder')
import VSM

def volume_var_mctigue(dpmu,r,d):
    r3 = r**3
    rd4 = (r/d)**4
    dVolume = np.pi*r3*dpmu*(1+rd4)
    return dVolume

def volume_var_yang(dpmu,s_axis_max,ratio):
    if(ratio > 1): ratio = 1./ratio
    volume = 4./3.*np.pi*s_axis_max**3*ratio*ratio
    f = ratio*ratio/3. - 0.7*ratio + 1.37
    dVolume = volume*dpmu*3./4.*f
    print('factor',f)
    return dVolume,f

def volume_var_penny(dpmu,r,nu):
    r3 = r**3
    dVolume = 8./3.*(1-nu)*r3*dpmu
    return dVolume

def fault_vertices(TLC,l,w,strike,dip):
    s  = np.deg2rad(strike)
    cs = np.cos(s)
    ss = np.sin(s)
    d  = np.deg2rad(dip)
    cd = np.cos(d)
    sd = np.sin(d)
    vertices = np.zeros((4,3))
    vertices[0] = TLC
    trace = np.zeros((2,3))
    
    # Bottom Left Corner
    coox = TLC[0] + w*cd*cs
    cooy = TLC[1] - w*cd*ss
    cooz = TLC[2] + w*sd
    vertices[1] = [coox,cooy,cooz]
	
    # Bottom Right Corner
    coox = vertices[1,0] + l*ss
    cooy = vertices[1,1] + l*cs
    cooz = vertices[1,2]
    vertices[2] = [coox,cooy,cooz]

    # Up Right Corner
    coox = TLC[0] + l*ss
    cooy = TLC[1] + l*cs
    cooz = TLC[2]
    vertices[3] = [coox,cooy,cooz]

    # Trace Left
    ww = TLC[2]/sd
    coox = TLC[0] - ww*cd*cs
    cooy = TLC[1] + ww*cd*ss
    cooz = 0.0
    trace[0] = [coox,cooy,cooz]
    
    # Trace Right
    ww = TLC[2]/sd
    coox = vertices[3,0] - ww*cd*cs
    cooy = vertices[3,1] + ww*cd*ss
    cooz = 0.0
    trace[1] = [coox,cooy,cooz]

    return vertices,trace
  
if(__name__ == "__main__"):

    
#----------------------------------------------------------------------
# FAULT or DIKE VERTICES from the TOP LEFT CORNER (TLC)
    '''
    TLC = [7.48873e+05,9.82e+06,2.21526e+03]
    length = 9000.
    width  = 4000.
    strike = 180.
    dip    = 83.
    
    vert,trace = fault_vertices(TLC,length,width,strike,dip)
    print('Fault vertices ordered as TLC - BLC - BRC - TRC',vert)
    print('Trace from left to right', trace)
    '''
    
    
#----------------------------------------------------------------------
# VOLUME VARIATION SPHERE (MCTIGUE)
    '''
    a = 600
    dP_mu = 7.47e-3
    depth = 2612
    
    dVol = volume_var_mctigue(dP_mu,a,depth)
    print('Volume variation of the finite volume sphere',dVol)
    '''
    
    
#----------------------------------------------------------------------
# VOLUME VARIATION SPHEROID (YANG)
    '''
    a = 8600
    dP_mu = 1e-3
    ratio = 0.33
    dVol,f = volume_var_yang(dP_mu,a,ratio)
    print('Volume variation of the finite volume spheroid',dVol/1e6, '10e6 m3')
    '''

#----------------------------------------------------------------------
# VOLUME VARIATION SILL (FIALKO) in point-source approximation
# holds for depth/radius > 2 or 2.5
    '''
    r = 1000
    dP_mu = 3.34987e-3
    dVol = volume_var_penny(dP_mu,r,0.25)
    print('Volume variation of the penny-shaped crack',dVol/1e6,' 10e6 m3')
    '''
    
    
#----------------------------------------------------------------------
# PLOT 1D - 2D statistics via corner
    
    folder = '../USECASE/mogi1NA/'
    f_models ='VSM_models.csv'
    f_truths = 'VSM_best.csv'
    f_plot = 'VSM_plot1D2D_new.pdf'
    nskip = 4000
    
    
    df = pd.read_csv(os.path.join(folder,f_models))
    tf = pd.read_csv(os.path.join(folder,f_truths))
    nunk = np.size(tf)
    
    if(f_truths[4:8] == 'best'):
        dd = df.iloc[:,2:] # NA
        tt = tf.values.reshape(nunk) #NA
    else:
        dd = df #BI
        tp = tf.values.reshape(nunk)
        tt = tp[::2] #BI
    
    VSM.plot_1D_2D(dd, tt, dd.columns, os.path.join(folder,f_plot), nskip)
    

    
#----------------------------------------------------------------------
# Plot Data, models and residuals
'''
folder_inout = '../USECASE/mogi1NA'
synth_file = 'VSM_synth_sar2.csv'
out_file = 'VSM_res_sar2.png'

# UTM zone
zone = "33N"
southern_hemisphere = False

# My data in UTM coordinates
db_sar = pd.read_csv(os.path.join(folder_inout,synth_file))
d_sar = db_sar.values

# Split into east and north coordinates
east, north = d_sar[:,0],d_sar[:,1]
data = d_sar[:,3]
synth = d_sar[:,2]
res = data - synth

dmax = max(max(data),max(synth))
dmin = min(min(data),min(synth))
resmax = max(res)
resmin = min(res)
resext = max(resmax, -resmin)

# Define the projection
#crs=ccrs.PlateCarree()
mycrs=ccrs.UTM(zone=zone, southern_hemisphere=southern_hemisphere)

fig=plt.figure(figsize=(15,6))
palette='viridis' #'RdYlBu_r'

## PANEL DATA ##########
ax = plt.subplot(131, projection=mycrs)
ax.coastlines(resolution='10m')
img = ax.scatter(east, north,5, data,cmap=palette, vmin=dmin, vmax = dmax)
#palette
cbar = plt.colorbar(img,orientation='horizontal')
cbar.set_label('LOS (m)')
# Get the extent of the axis
extent = ax.get_extent()
# Attempt to set the axis extent
ax.set_extent(extent, crs=mycrs)
plt.title('Data',fontsize = 16, pad=10)

## PANEL MODEL ##########
ax = plt.subplot(132, projection=mycrs)
ax.coastlines(resolution='10m')
img = ax.scatter(east, north,5, synth,cmap=palette, vmin=dmin, vmax = dmax)
#palette
cbar = plt.colorbar(img,orientation='horizontal')
cbar.set_label('LOS (m)')
# Get the extent of the axis
extent = ax.get_extent()
# Attempt to set the axis extent
ax.set_extent(extent, crs=mycrs)
plt.title('Model',fontsize = 16, pad=10)

## PANEL RESIDUALS ##########
ax = plt.subplot(133, projection=mycrs)
ax.coastlines(resolution='10m')
img = ax.scatter(east, north,5, res,cmap="bwr",vmin=resext, vmax = -resext)
#palette
cbar = plt.colorbar(img,orientation='horizontal')
cbar.set_label('LOS (m)')
# Get the extent of the axis
extent = ax.get_extent()
# Attempt to set the axis extent
ax.set_extent(extent, crs=mycrs)
# Title for plot
plt.title('Residual',fontsize = 16, pad=10)

plt.savefig(os.path.join(folder_inout,out_file))
plt.show()
'''