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
"""

#import sys
#sys.path.append('path-to-VSM-folder')
import VSM

import tkinter as tk
from   tkinter import filedialog
from   tkinter import ttk

fontstring = 'Helvetica 13'
modcolor = 'lightgrey'

def letsgo():
    thatsallfolks(r)
    
    form = tk.Tk()
    form.title('VSM Input')
    form.geometry('470x460')
    
    note = ttk.Notebook(form)
    
    tab1= ttk.Frame(note)
    tab2= ttk.Frame(note)
    tab3= ttk.Frame(note)
    
    note.add(tab1, text = 'Data')
    note.add(tab2, text = 'Model')
    note.add(tab3, text = 'Inversion')
    note.pack(expand=1,fill='both')
    
    #TAB1
    GeoData = ttk.LabelFrame(tab1, text = "Geodetic Data Files & Weights")
    GeoData.grid(column=0, row=0, sticky='E', padx=8, pady=4) 
    
    ButtonRead = ttk.Button(GeoData, text = 'InSAR Input File(s)', command = VSM.read_sar_file)
    ButtonRead.grid(column = 0, row = 1, sticky='W',padx = 8,pady=5)
    
    ButtonRead = ttk.Button(GeoData, text = 'GNSS Input File', command = VSM.read_gps_file)
    ButtonRead.grid(column = 0, row = 2, sticky='W',padx = 8,pady=5)
    
    ButtonRead = ttk.Button(GeoData, text = 'Levelling Input File', command = VSM.read_lev_file)
    ButtonRead.grid(column = 0, row = 3, sticky='W', padx = 8,pady=5)
    
    ButtonRead = ttk.Button(GeoData, text = 'EDM Input File', command = VSM.read_edm_file)
    ButtonRead.grid(column = 0, row = 4, sticky='W', padx = 8,pady=5)
    
    ButtonRead = ttk.Button(GeoData, text = 'Tilt Input File', command = VSM.read_tlt_file)
    ButtonRead.grid(column = 0, row = 5, sticky='W', padx = 8,pady=5)

    ButtonRead = ttk.Button(GeoData, text = 'Strain Input File', command = VSM.read_srn_file)
    ButtonRead.grid(column = 0, row = 6, sticky='W', padx = 8,pady=5)

    ttk.Label(GeoData, text = "InSAR Weight").grid(    column=2,row=1,sticky='E')
    ttk.Label(GeoData, text = "GNSS Weight").grid(     column=2,row=2,sticky='E')
    ttk.Label(GeoData, text = "Levelling Weight").grid(column=2,row=3,sticky='E')
    ttk.Label(GeoData, text = "EDM Weight").grid(      column=2,row=4,sticky='E')
    ttk.Label(GeoData, text = "Tilt Weight").grid(     column=2,row=5,sticky='E')
    ttk.Label(GeoData, text = "Strain Weight").grid(   column=2,row=6,sticky='E')
    
    WG_sar = tk.Entry(GeoData, width = 6, bg = 'Lightblue')
    WG_sar.insert(0,"0")
    WG_sar.focus_set()
    WG_sar.grid(column = 3, row = 1)
    WG_gps = tk.Entry(GeoData, width = 6, bg = 'Lightblue')
    WG_gps.insert(0,"0")
    WG_gps.focus_set()
    WG_gps.grid(column = 3, row = 2)
    WG_lev = tk.Entry(GeoData, width = 6, bg = 'Lightblue')
    WG_lev.insert(0,"0")
    WG_lev.focus_set()
    WG_lev.grid(column = 3, row = 3)
    WG_edm = tk.Entry(GeoData, width = 6, bg = 'Lightblue')
    WG_edm.insert(0,"0")
    WG_edm.focus_set()
    WG_edm.grid(column = 3, row = 4)
    WG_tlt = tk.Entry(GeoData, width = 6, bg = 'Lightblue')
    WG_tlt.insert(0,"0")
    WG_tlt.focus_set()
    WG_tlt.grid(column = 3, row = 5)
    WG_srn = tk.Entry(GeoData, width = 6, bg = 'Lightblue')
    WG_srn.insert(0,"0")
    WG_srn.focus_set()
    WG_srn.grid(column = 3, row = 6)
    
    for i in range(5): ttk.Label(tab1,text = '').grid()
    tk.Button(tab1, text='EXIT', command = lambda : (thatsallfolks(form),exit())).grid()             

    #TAB2  
    ModData = ttk.LabelFrame(tab2, text = "Elastic Constants")
    ModData.grid(padx=8, pady=4,column = 0, sticky='W')

    ttk.Label(ModData, text = "Shear Modulus (Pa)").grid(column=0,row=1,padx=8, pady=4,sticky='W')
    ttk.Label(ModData, text = "Poisson Coefficient").grid(column=0,row=2,padx=8, pady=4,sticky='W')

    MM = tk.Entry(ModData, width = 6, bg = 'Lightblue')
    MM.insert(0,"5e9")
    MM.focus_set()
    MM.grid(column = 1, row = 1,sticky = 'W')
    NN = tk.Entry(ModData, width = 6, bg = 'Lightblue')
    NN.insert(0,"0.25")
    NN.focus_set()
    NN.grid(column = 1, row = 2,sticky = 'W')
    
    ttk.Label(ModData, text = "Forward Model").grid(column=0,row=3,padx=8, pady=4,sticky='W')

    txt=[]
    for i in range(6): txt.append(VSM.source_info(i)['name'])
    
    ButtonRead = tk.Button(ModData, text = txt[0], 
                           command = lambda: read_model(0))
    ButtonRead.grid(row = 4, sticky='W',padx = 8,pady=5)
    
    ButtonRead = tk.Button(ModData, text = txt[1], 
                           command = lambda: read_model(1))
    ButtonRead.grid(row = 5, sticky='W',padx = 8,pady=5)
    
    ButtonRead = tk.Button(ModData, text = txt[2], 
                           command = lambda: read_model(2))
    ButtonRead.grid(row = 6, sticky='W',padx = 8,pady=5)
    
    ButtonRead = tk.Button(ModData, text = txt[3], 
                           command = lambda: read_model(3))
    ButtonRead.grid(row = 7, sticky='W',padx = 8,pady=5)
    
    ButtonRead = tk.Button(ModData, text = txt[4], 
                           command = lambda: read_model(4))
    ButtonRead.grid(row = 8, sticky='W',padx = 8,pady=5)
    
    ButtonRead = tk.Button(ModData, text = txt[5], 
                           command = lambda: read_model(5))
    ButtonRead.grid(row = 9, sticky='W',padx = 8,pady=5)
    
    ttk.Label(tab2,text = '').grid(row=10)
    tk.Button(tab2, text='EXIT',
          command = lambda : (thatsallfolks(form),exit())).grid(row=11, column=0, pady=4)             
    
    #TAB3
    OPTData = tk.LabelFrame(tab3, text = "Inversion Method",font=fontstring)
    OPTData.grid(column=0, row=1, padx=4, pady=4) 

    OptionList = ["Neighbourhood Algorithm", "Bayesian Inference"]
    variable = tk.StringVar(OPTData)
    variable.set(OptionList[0])

    opt0 = tk.Radiobutton(OPTData,text=OptionList[0], value=0,command = lambda: read_opt(0))
    opt0.config(width=45, font=fontstring)
    opt0.grid(column = 0, row = 1, sticky = 'W')
    opt1 = tk.Radiobutton(OPTData,text=OptionList[1], value=1,command = lambda: read_opt(1))
    opt1.config(width=45, font=fontstring)
    opt1.grid(column = 0, row = 2, sticky = 'W')
        

    OUTData = ttk.LabelFrame(tab3, text = "Output Folder")
    OUTData.grid(column=0, row=2, sticky='W', padx=8, pady=4) 

    ButtonRead = tk.Button(OUTData, text = 'In/Out Folder', command = lambda: VSM.read_inout_folder())
    ButtonRead.grid(column = 1, row = 2, sticky='W',padx = 8,pady=5)

    plotyn = ttk.LabelFrame(tab3, text = "Plots ")
    plotyn.grid(column=0, row=6, sticky='W',padx=8, pady=4)

    ttk.Label(plotyn, text = "Number of Burn-in Models").grid(column=0,row=7,sticky='W')

    BI = tk.Entry(plotyn, width = 6, bg = 'Lightblue')
    BI.insert(0,"2000")
    BI.focus_set()
    BI.grid(column = 1, row = 7)
     
    ttk.Label(tab3,text = '').grid()
    tk.Button(tab3, text='Launch VSM!',
          command = lambda : (VSM.read_weights(WG_sar.get(),WG_gps.get(),WG_lev.get(),WG_edm.get(),WG_tlt.get(),WG_srn.get()),
                              VSM.read_elast(MM.get(),NN.get()), VSM.read_plot(BI.get()),
                              VSM.final_check(),thatsallfolks(form))).grid()
    
    for i in range(4): ttk.Label(tab3,text = '').grid()
    tk.Button(tab3, text='EXIT',
          command = lambda : (thatsallfolks(form),exit())).grid()                                                                                     
    
    form.mainloop()
    
def preread_optim(option,entries):
    
    mm = []
    for entry in entries:
        mm.append(entry.get())
    
    if(option == 1): mm.append(10)
        
    VSM.read_optim(option,*mm)

def preread_params(entries1,entries2,opt):
    m1 = []
    m2 = []
    for entry in entries1:
        m1.append(entry.get())
    for entry in entries2:
        m2.append(entry.get())
    VSM.read_params(m1,m2,opt)
    
def read_lic():
    rlic = tk.Tk()
    LicArea = tk.Frame(rlic, bg = modcolor)
    LicArea.grid()
    tk.Label(LicArea, text='Disclaimer',bg = modcolor).grid()
    textbox = tk.Text(LicArea) 
    flic = open("./license.lic", "r")
    i = 1.0
    for x in flic:
        textbox.insert(i, x+'\n')
        i+=1
    textbox.grid()
    tk.Button(LicArea, text = 'Done', command = lambda: thatsallfolks(rlic)).grid()
    rlic.mainloop() 

def read_model(option):
    global model
    
    model = VSM.source_info(option)

    rmod = tk.Tk()
    txt = model['name']
    rmod.title(txt)
    rmod.config(bg = modcolor)
 
 
    shift = 0
    
    if(option == 5):
        OPTMod = tk.LabelFrame(rmod, text = "Select Okada mode",font=fontstring,bg=modcolor)
        OPTMod.grid(column=0, row=1, padx=6, pady=4) 
    
        OptionList = ["Total slip (param1) and rake (param2)", "Strike slip (param1) and dip slip (param2)"]
        variable = tk.StringVar(OPTMod)
        variable.set(OptionList[0])

        opt0 = tk.Radiobutton(OPTMod,text=OptionList[0], value=0,command = lambda: VSM.read_Okada_mode('R'))
        opt0.config(width=45, font=fontstring,bg=modcolor)
        opt0.grid(column = 0, row = 1, sticky='W')
        opt1 = tk.Radiobutton(OPTMod,text=OptionList[1], value=1,command = lambda: VSM.read_Okada_mode('S'))
        opt1.config(width=45, font=fontstring,bg=modcolor)
        opt1.grid(column = 0, row = 2, sticky='W')
        shift = 3

    ParArea = tk.Frame(rmod,bg = modcolor)
    ParArea.grid()
    
    tk.Label(ParArea, text = "Parameter",font = fontstring,bg=modcolor).grid(column=1,row=shift)
    tk.Label(ParArea, text = "Min",      font = fontstring,bg=modcolor).grid(column=2,row=shift)
    tk.Label(ParArea, text = "Max",      font = fontstring,bg=modcolor).grid(column=3,row=shift)
        
    en1 = []
    en2 = []
    for i in range(model['param_n']):
        
        tk.Label(ParArea, text = model['param_l'][i],bg=modcolor,font =fontstring).grid(column=1,row=shift+i+1)
        M1 = tk.Entry(ParArea, width = 14, bg = 'Lightblue')
        M1.insert(0,"")
        M1.focus_set()
        M1.grid(column = 2, row = shift+i+1,sticky = 'W')
        en1.append(M1)
        M2 = tk.Entry(ParArea, width = 14, bg = 'Lightblue')
        M2.insert(0,"")
        M2.focus_set()
        M2.grid(column = 3, row = shift+i+1,sticky = 'W')
        en2.append(M2)
        tk.Button(ParArea, text='Done!',
             command = lambda : (preread_params(en1,en2,option),thatsallfolks(rmod))).grid(
                             row=shift+i+2, column=2, pady=4,sticky='W')

    rmod.mainloop()

def read_opt(option):
    rmod = tk.Tk()
    txt = 'Choose set up'
    rmod.title(txt)
    rmod.config(bg = modcolor)

    
    en = []
    if(option == 0):
        ParArea = tk.LabelFrame(rmod, text = "Neighbourhood Algorithm parameters", bg = modcolor)
        ParArea.grid(column=0, row = 0)
        txt = ['N. of samples','N. of resamples','N. of iterations']
        init = ["1000","300","10"]
        for i in range(3):
            tk.Label(ParArea, text = txt[i],bg=modcolor).grid(column=0,row=i+1,sticky='W',padx=3)
            M1 = tk.Entry(ParArea, width = 8, bg = 'Lightblue')
            M1.insert(0,init[i])
            M1.focus_set()
            M1.grid(column = 1, row = i+1,sticky = 'W')
            en.append(M1)
    if(option == 1):
        ParArea = tk.LabelFrame(rmod, text = "Bayesian Inference parameters", bg = modcolor)
        ParArea.grid(column=0, row = 0)
        txt = ['N. of random walks','N. of steps']
        init = ["18","8000"]
        for i in range(2):
            tk.Label(ParArea, text = txt[i],bg=modcolor).grid(column=0,row=i+1,sticky='W',padx=3)
            M1 = tk.Entry(ParArea, width = 8, bg = 'Lightblue')
            M1.insert(0,init[i])
            M1.focus_set()
            M1.grid(column = 1, row = i+1,sticky = 'W')
            en.append(M1)
    
    tk.Button(ParArea, text='Done!',
             command = lambda : (preread_optim(option,en),thatsallfolks(rmod))).grid(
                             row=i+2, column=0, pady=4,sticky='W')

    rmod.mainloop()
    
def read_VSM_settings_all():
    
    VSM_settings_file = filedialog.askopenfilename(title = "Select File",filetypes = (("ascii files","*.txt"),("all files","*.*")))
    VSM.read_VSM_settings(VSM_settings_file)
    thatsallfolks(r)

def thatsallfolks(name):
    name.destroy()
    
if __name__ == "__main__":

    vsmcolor='lightgrey'
    r = tk.Tk()
    r.title('VSM - Volcanic and Seismic source Modelling')
    r.config(bg = vsmcolor)
    
    logo = tk.PhotoImage(file ="./VSM_logo.gif")
    logo_label = tk.Label(image = logo)
    logo_label.pack(padx=1)
    
    VSMInit = tk.Label(r, text = "", font = ('Helvetica 16'), bg = vsmcolor)
    VSMInit.pack(padx = 5)


    ButtonArea = tk.Frame(r, bg = vsmcolor)
    ButtonArea.pack()
    
    ButtonRead = tk.Button(ButtonArea, text = 'Load VSM Input', bg = vsmcolor, command = read_VSM_settings_all)
    ButtonRead.grid(column = 0, row = 0, padx = 5)
    
    
    ButtonInsert = tk.Button(ButtonArea, text = 'Create VSM Input', bg = vsmcolor, command = letsgo)
    ButtonInsert.grid(column = 1, row = 0, padx = 5)


    ButtonInsert = tk.Button(ButtonArea, text = 'License', bg = vsmcolor, command = read_lic)
    ButtonInsert.grid(column = 2, row = 0, padx = 5)

    
    ButtonExit = tk.Button(ButtonArea, text = 'EXIT', bg = vsmcolor, command = lambda: (thatsallfolks(r), exit()))
    ButtonExit.grid(column = 3, row = 0, padx = 5,pady=4)
    
    BelowArea = tk.Label(r, text = "", bg = vsmcolor)
    BelowArea.pack(pady=5)
        
    
    r.mainloop()
    
    VSM.iVSM()