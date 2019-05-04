
# coding: utf-8

# In[ ]:


import numpy as np
import os


# In[ ]:


def MAP_mat(Directorylocation):  #returns matrix with shape (xspot number, yspot number )
    nXspot=1 
    nYspot=1
    stepX=0
    stepY=0
    nPixels=2048
    #Directorylocation='Hydro Al sheets_120 pulse_10pbc_1ps/5Z-0_1x10x120/'
    txtNames=os.listdir(Directorylocation)
    #######################################
    #cheching the  number of pulses, steps and spots
    for i in range(len(txtNames)-1):
        if int(txtNames[i][6:12])!=int(txtNames[i+1][6:12]) or int(txtNames[i][14:20])!=int(txtNames[i+1][14:20]):
            nPulsePerSpot=i+1
            break     
    for i in range(len(txtNames)-1):
        if int(txtNames[i][6:12])!=int(txtNames[i+1][6:12]):
            stepY=int(txtNames[i+1][6:12])-int(txtNames[i][6:12])
            nYspot+=1
        if int(txtNames[i][14:20])!=int(txtNames[i+1][14:20]):
            stepX=int(txtNames[i+1][14:20])-int(txtNames[i][14:20])
            nXspot+=1
    ################################################
    #initialization
    yInitialPosition=int(txtNames[0][6:12])
    xInitialPosition=int(txtNames[0][14:20])
    c=np.asarray([txtNames[i][22:len(txtNames[i])-4] for i in range(len(txtNames))])
    Ccount=0
    s='000000'
    info=np.zeros((nXspot,nYspot,nPulsePerSpot,nPixels))
    ####################################
    #Construction
    
    for ny in range(nYspot):
        for nx in range(nXspot):
            #updating positions 
            yposition=yInitialPosition+ny*stepY
            xposition=xInitialPosition+nx*stepX
            
            for p in range(nPulsePerSpot):
                #constructing file name to be read
                Ynumber=s[:6-len(str(yposition))]+str(yposition)
                Xnumber=s[:6-len(str(xposition))]+str(xposition)
                
                filename='spec_Y'+Ynumber+'_X'+ Xnumber+'_C'+c[Ccount]+'.txt'
                
                DataTXT=np.loadtxt(Directorylocation+filename, delimiter="\t",skiprows=1)
                info[nx,ny,p]=DataTXT[:,3]
                Ccount+=1
                
        
           
    return info

