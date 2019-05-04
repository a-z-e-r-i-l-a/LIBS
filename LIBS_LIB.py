
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from matplotlib import pyplot, pylab
import numpy as np
import os
import seaborn as sns
import pandas as pd
from scipy import signal
from scipy.misc import electrocardiogram
from scipy.signal import find_peaks
from scipy import sparse
from scipy.sparse.linalg import spsolve


# In[2]:


def MAP_mat(Directorylocation,bg_removal=0):  #returns matrix with shape (xspot number, yspot number )
    
    nXspot=1 
    nYspot=1
    stepX=0
    stepY=0
    nPixels=2048
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
    #
    nXspot=1#int(nXspot/nYspot)
    ################################################
    #initialization
    yInitialPosition=int(txtNames[0][6:12])
    xInitialPosition=int(txtNames[0][14:20])
    c=np.asarray([txtNames[i][22:len(txtNames[i])-4] for i in range(len(txtNames))])
    Ccount=0
    s='000000'
    ####################################
    #Construction
    info=np.zeros((nXspot,nYspot,nPulsePerSpot,nPixels))

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

    #wavelengths=DataTXT[:,0]

    if bg_removal==1:
        spectra=[]
        for i in range(info.shape[2]):
            spectra.append(baseline_removal(np.squeeze(spectra)[:,i,:]))
    else:
        spectra=info


    return np.asarray(spectra)

def t_resolved_spectra(y_coordinate,x_coordinate,Directorylocation,nPulsePerSpot=100,nPixels=2048):
  

    txtNames=os.listdir(Directorylocation)
    c=np.asarray([txtNames[i][22:len(txtNames[i])-4] for i in range(len(txtNames))])
    all_spots_coordinates=np.asarray([txtNames[i][5:20] for i in range(len(txtNames))])

    info=np.zeros((nPulsePerSpot,nPixels))

    spot_coordinate='Y'+'000000'[:6-len(str(y_coordinate))]+str(y_coordinate)+'_X'+'000000'[:6-len(str(x_coordinate))]+str(x_coordinate)
    #current_coordinate='Y'+'000000'[:6-len(str(y_coordinate))]+str(y_coordinate)+'_X'+'000000'[:6-len(str(x_coordinate))]+str(x_coordinate)
    C_idx=np.where(all_spots_coordinates==spot_coordinate)

    for p in range(nPulsePerSpot):
        #constructing file name to be read

        filename='spec_'+spot_coordinate+'_C'+c[C_idx][p]+'.txt'

        DataTXT=np.loadtxt(Directorylocation+'/'+filename, delimiter="\t",skiprows=1)
        info[p]=DataTXT[:,3]


    wavelengths=DataTXT[:,0]
    return wavelengths, info
###############################################################################################################

def AllTXTfilesSpectrumMatrix(location):
    os.chdir(location)
    
    txtNames=[]
    for file in os.listdir(location):
        if file.endswith(".txt"):
            if np.loadtxt(file, delimiter="\t",skiprows=0,dtype='str')[0,3]=='SpecMinusDark': #making sure its spectroscope txt file
                txtNames.append(file)

    intensities=np.zeros((len(txtNames),2048))
    for i in range(len(txtNames)):
        #file = open(txtNames[i], encoding="utf8")
        intensities[i,:]=np.loadtxt(txtNames[i], delimiter="\t",skiprows=1,encoding='bytes')[:,3]
    wavelengths=np.loadtxt(txtNames[0], delimiter="\t",skiprows=1)[:,0] 
    return txtNames, wavelengths, intensities
# give locaion of the folder full of txt files from the avantes spectrometer,
#it retruns the txt names in order
#it returns the wavelengths
#it returns a matrix of Nx2048 from all txt files

############################################################################################################################

def baseline(y, lam=0.1, p=0.001, niter=10): #1x2048
    
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z #1x2048

def baseline_removal(spectra,lam=0.1, positive=0.001, niteration=10): #Nx2048
        
    spectra=np.atleast_2d(spectra)
    spectra_without_BL=np.zeros_like(spectra)
    for i in range(len(spectra)):
        spectra_without_BL[i]=spectra[i]-baseline(spectra[i],lam, positive, niteration)
    return np.squeeze(spectra_without_BL)  #Nx2048
############################################################################################################################

def draw_heatmap(spectra_without_BL,wavelengths,names,label='Spectra without the baseline'): #Nx2048

    spectra_without_BL=pd.DataFrame(spectra_without_BL)
    spectra_without_BL.columns=wavelengths
    spectra_without_BL.index=names

    plt.figure(figsize=[50,20])
    sns.heatmap(spectra_without_BL,xticklabels = 25, yticklabels = 1 ,robust=1, cbar_kws={'label': label})
    
def give_peaks_loc(spectra,wavelength,p_height=1000,peakdistance=1,lam=0.1, positive=0.001, niterations=10,disp=1):
    spectra_no_BG=baseline_removal(spectra,lam, positive, niterations)
    from scipy.signal import find_peaks
    
    loc=find_peaks(spectra_no_BG,distance=peakdistance,height=p_height)[0][:]
    
    if disp==1:
        plt.figure(figsize=[25,5])
        f=plt.plot(wavelengths,spectra)
        f=plt.scatter(wavelengths[loc],spectra[loc],color='red')

        plt.figure(figsize=[25,5])
        f=plt.plot(wavelengths,spectra_no_BG)
        f=plt.scatter(wavelengths[loc],spectra_no_BG[loc])
        
        
    return loc #location of the peaks in the wavelength vector
    
    
def find_best_non_saturated_pixels(spectra,hight_threshold=200):
    Background_removed_spectra=baseline_removal(spectra)
    w=np.ones_like(wavelengths)
    for n in range(len(spectra)):
        #inclusion criteria
        idx_with_1st_crit=np.where(spectra[n]<16000)[0]
        idx_with_2nd_crit=np.where(Background_removed_spectra[n]>=hight_threshold)[0]
        exclusionCrit=np.intersect1d(idx_with_1st_crit,idx_with_2nd_crit)
        #exclusion criterion
        exclusion_idxs=[i for i in range(len(w)) if i not in exclusionCrit]
        w[exclusion_idxs]*=0
        idx_with_1st_crit=None
        idx_with_2nd_crit=None
        exclusion_idxs=None
    plt.figure(figsize=[30,5])
    plt.plot(wavelengths,w)
    
    loc=np.where(w==1)[0]
    Wavelengths=wavelengths[loc]
    std=np.sqrt(np.var(Background_removed_spectra[:,loc],axis=0))
    mean=np.mean(Background_removed_spectra[:,loc],axis=0)  
    minv=np.min(Background_removed_spectra[:,loc],axis=0)  
    maxv=np.max(Background_removed_spectra[:,loc],axis=0)  
    
    info=pd.DataFrame(np.vstack((Wavelengths,loc,mean,minv,maxv,std)))
    info.index='Wavelength','Pixel number','Mean value','Minimum value','Maximum value','Standard deviation'
    display(info)
    return loc, Wavelengths, std
#gives the pixels that are not saturated in any of the 
#spectra and are higher than the threshold in the baseline-line removed of all spectra



def find_nearest(array, value):
    value=np.atleast_1d(np.asarray(value))
    array = np.asarray(array)
    idx=[]
    for i,value_i in enumerate(value):
        loc=(np.abs(array - value_i)).argmin()
        if loc not in idx:
            idx.append(loc)
    return idx,array[idx]

def spectra2new_framework(spectra,initial_wavelength=200,final_wavelength=800,resolution=0.5,p_distance=1,peak_height=2000):
    #removing the baseline
    spectra_no_BG=baseline_removal(spectra)
    ####################################################################################################
    wavelengths_range=np.arange(initial_wavelength,final_wavelength+resolution,resolution)
    peaks_values_in_new_frame=np.zeros((len(spectra_no_BG),len(wavelengths_range)))
    for s in range(len(spectra_no_BG)):
        p_idx=find_peaks(spectra_no_BG[s],height=peak_height,distance=p_distance)[0][:]
        for p,p_loc in enumerate(p_idx):
            near_loc,_=find_nearest(wavelengths_range,wavelengths[p_loc])
            peaks_values_in_new_frame[s,near_loc]=spectra_no_BG[s,p_loc]

        P_loc=None
    return spectra_no_BG,peaks_values_in_new_frame

def X2Fmap(X): ## X array with Dx1 dimansion
    idx=np.where(X!=0)[0]
    return idx.reshape(-1,1)


# In[3]:


def spectra_peaks(y):
    spectra_peaks=np.zeros_like(y)
    for i in range(len(y)):
        peak_idx =  find_LIBS_peaks  (y[i],peak_distance=1,peak_height =1000,threshold=1,closeness=10,plot=0)
        spectra_peaks[i,peak_idx]=1
        peak_idx=None
    return spectra_peaks# similar to y only the peaks with value 1


# In[4]:


def baseline(y, lam=0.1, p=0.001, niter=10): #1x2048
    
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L), 2))
    w = np.ones(L)
    for i in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z #1x2048

def baseline_removal(spectra,lam=0.1, positive=0.001, niteration=10): #Nx2048
        
    spectra=np.atleast_2d(spectra)
    spectra_without_BL=np.zeros_like(spectra)
    for i in range(len(spectra)):
        spectra_without_BL[i]=spectra[i]-baseline(spectra[i],lam, positive, niteration)
    return np.squeeze(spectra_without_BL)  #Nx2048


# In[5]:


def Normalized_spectra(spectra):
    spectra_normalized=np.zeros_like(spectra)
    for i,spectrum in enumerate(spectra):
        spectra_normalized[i]=spectrum/np.max(spectrum)
    return spectra_normalized

def unit_peaks_spectra(spectra,normalize=1,peak_threshold=0.25):
    if normalize==1:
        spectra1=Normalized_spectra(spectra)
    spectra2=np.zeros_like(spectra)
    for i in range(len(spectra)):
        idx=np.where(spectra[i]>=peak_threshold)
        spectra2[i,idx]=1
    return spectra2


# In[6]:


def create_sorted_training_set(location,peak_threshold=1000,lam=0.001, positive=0.000001, niteration=20):

    txtNames, wavelengths, spectra=AllTXTfilesSpectrumMatrix(location)
    spectra_no_BG=baseline_removal(spectra,lam=lam, positive=positive, niteration=niteration)
    file_names=np.asarray([name[:2] for i,name in enumerate(txtNames)])
    ###############################   SORTING THE SPECTRA BASED ON THE NUMBER OF PEAKS
    NumberOfPeaks=np.zeros(len(file_names))
    for i in range(len(file_names)):
        idx=find_peaks(spectra_no_BG[i],distance=1,height=peak_threshold)[0][:]
        NumberOfPeaks[i]=idx.shape[0]

    sorted_idx=NumberOfPeaks.argsort()[::-1] #indecies in decending order

    NumberOfPeaks=NumberOfPeaks[sorted_idx]
    txtNames_sorted=np.asarray(txtNames)[sorted_idx]
    spectra_sorted=spectra[sorted_idx]
    spectra_no_BG_sorted=spectra_no_BG[sorted_idx]
    spectra_w_BG_sorted=spectra[sorted_idx]
    file_names_sorted=file_names[sorted_idx]
    _, idx=np.unique(file_names_sorted,return_index=True)
    Material_names_sorted=file_names_sorted[np.sort(idx)]
    ###   assinging a one-hot vector to each material calsss    ######
    n_file=len(file_names_sorted)
    n_class=len(Material_names_sorted)
    One_hot_vector=np.zeros((1,len(Material_names_sorted)))

    y_train=np.zeros((n_file,n_class))
    for i in range(n_class):
        One_hot_vector=np.zeros((1,n_class))
        One_hot_vector[0,i]=1
        idx=np.where(file_names_sorted==Material_names_sorted[i])[0]
        y_train[idx,:]=One_hot_vector

    X_train=spectra_no_BG_sorted

    return X_train, y_train, file_names_sorted, Material_names_sorted,wavelengths
#X_train 54x2048 
#y_train 54x15


# In[7]:


def Marked_peaks_plot(spectrum,wavelengths,peak_threshold=150,name=None,single_plot=1):

    idx=find_peaks(spectrum,distance=2,height=peak_threshold)[0][:]
    
    if len(idx) is not 0:
        idx_sorted=idx[spectrum[idx].argsort()][::-1]
        sorted_peaks=pd.DataFrame(np.column_stack((np.atleast_2d(wavelengths[idx_sorted]).T,
                                                   np.atleast_2d(spectrum[idx_sorted]).T,
                                                   np.atleast_2d(idx_sorted).T)))
        sorted_peaks.columns=['Wavelength','Intensity','Pixel']


    ######## PLotting

        df=pd.DataFrame(w)

        #fig, ax = plt.subplots()
        if single_plot==1:
            plt.figure(figsize=[30,10])

        plt.plot(wavelengths,spectrum,linewidth=1,label='Spectrum without baseline')
        plt.scatter(wavelengths[idx],spectrum[idx],linewidth=1,label='Peaks',marker="v",color='m')
        ticks=plt.xticks(np.arange(min(wavelengths), max(wavelengths)+1, 10.0))

        for i, pixel in enumerate(idx_sorted):
            plt.annotate(str(i+1), (wavelengths[pixel]-0.8, spectrum[pixel]*1.025),
                         size=9)

        #### table
        
        
        
        row=[]
        w=[]
        valPerCol=45
        if len(idx_sorted)<valPerCol:
            numOfcolumns=1
        elif len(idx_sorted)>valPerCol and len(idx_sorted)<valPerCol*2:
            numOfcolumns=2
        elif len(idx_sorted)>valPerCol*2 and len(idx_sorted)<valPerCol*3:
            numOfcolumns=3
        elif len(idx_sorted)>valPerCol*3 and len(idx_sorted)<valPerCol*4:
            numOfcolumns=4
        else:
            numOfcolumns=5

        c=0
        for i,lam in enumerate(wavelengths[idx_sorted]):
            row.append(str(i+1)+' | '+str(lam)+' : '+str(idx_sorted[i]))
            c+=1
            if c==numOfcolumns:
                w.append(row)
                c=0
                row=[]
        df=pd.DataFrame(w)


        the_table=plt.table(cellText=df.values,loc='right',cellLoc='left')
        for key, cell in the_table.get_celld().items():
            cell.set_linewidth(0.3)
        the_table.set_fontsize(20)
        the_table.scale(0.0666*numOfcolumns, 1)
        
        return sorted_peaks
    else:
        error=name+' had no peak with the threshold '+str(peak_threshold)+'.'
        print(error)


# In[8]:


def make_plots(location,against_spectrum=None,peak_threshold=150,saving_folder_name='Plots'):
    subfolders=[]
    for root, dirs, files in os.walk(location):
        for name in files:
            if name.endswith((".txt")):
                subfolders.append(root)
                break

    for i in range(len(subfolders)):
        txtNames, wavelengths, spectra=AllTXTfilesSpectrumMatrix(subfolders[i])
        
        spectra_no_BG=baseline_removal(spectra,lam=0.001, positive=0.000001, niteration=20)
        if against_spectrum is not None:
            spectra_no_BG_2=baseline_removal(against_spectrum,lam=0.001, positive=0.000001, niteration=20)
            
        
        spectra=np.atleast_2d(spectra)
        spectra_no_BG=np.atleast_2d(spectra_no_BG)

        for i in range(len(txtNames)):
            plt.figure(figsize=[30,20])


            plt.subplot(2, 1,1)
            plt.plot(wavelengths,spectra[i],linewidth=1,label='Spectrum')
            if against_spectrum is not None:
                plt.plot(wavelengths,against_spectrum,linewidth=1,label='Refrence spectrum ')

            plt.xticks(np.arange(min(wavelengths), max(wavelengths)+1, 10.0))
            plt.title(txtNames[i])
            plt.legend()

            plt.subplot(2, 1,2)
            #plt.plot(wavelengths,spectra_no_BG[i],linewidth=1,label='Spectrum without baseline')
            sorted_peaks=Marked_peaks_plot(spectra_no_BG[i],wavelengths=wavelengths,peak_threshold=peak_threshold,name=txtNames[i],single_plot=0)
            
            if against_spectrum is not None:
                plt.plot(wavelengths,spectra_no_BG_2,linewidth=1,label='Refrence spectrum without baseline')

            plt.xticks(np.arange(min(wavelengths), max(wavelengths)+1, 10.0))
            plt.legend()

            
                #saving the plots
            if not os.path.exists(saving_folder_name):
                os.makedirs(saving_folder_name)
            plt.savefig(saving_folder_name+'/'+ txtNames[i]+'.jpg',dpi=250)
            plt.close()
            
            
                #Saving the peaks information
            if sorted_peaks is not None:
                saving_folder_name2='Peaks Information'
                if not os.path.exists(saving_folder_name2):
                    os.makedirs(saving_folder_name2)             
                sorted_peaks.to_csv(saving_folder_name2+'/'+ txtNames[i]+'.csv',header=True,index=False)
            sorted_peaks=None
            
            
    pass


# In[9]:


def find_libs_sorted_peaks(spectrum,wavelengths,peak_threshold=150):
    idx=find_peaks(spectrum,distance=2,height=peak_threshold)[0][:]
    
    if len(idx) is not 0:
        idx_sorted=idx[spectrum[idx].argsort()][::-1]
        sorted_peaks=pd.DataFrame(np.column_stack((np.atleast_2d(wavelengths[idx_sorted]).T,
                                                   np.atleast_2d(spectrum[idx_sorted]).T,
                                                   np.atleast_2d(idx_sorted).T)))
        sorted_peaks.columns=['Wavelength','Intensity','Pixel']
        
        return idx_sorted,sorted_peaks
    else:
        print('no peaks with the threshold '+str(peak_threshold))
        pass


# In[10]:


def Analyse_And_Read_Updated_Database(Avantes=1,
                                      peak_threshold=150,
                                      lam=0.001,
                                      positive=0.000001,
                                      niteration=20):
    if Avantes==1:
    
        location=r'\\ipmlan\PK\Projekte\Fraunhofer\600610_InlineElement\10_akt_Arbeiten_Mitarbeiter\Alireza\Automatic libs sample analysis with database\270_540nm'
        X_train, _, _, Material_names_sorted,wavelengths=create_sorted_training_set(location,
                                                                                      peak_threshold=peak_threshold,
                                                                                      lam=lam,
                                                                                      positive=positive,
                                                                                      niteration=niteration)
        ref_sorted_peaks=[]
        for i in range(len(X_train)):
            idx_sorted,sorted_peaks=find_libs_sorted_peaks(X_train[i],
                                                           wavelengths,
                                                           peak_threshold=peak_threshold)

            ref_sorted_peaks.append(sorted_peaks)
            ref_sorted_peaks[i]['Name']=Material_names_sorted[i]

        ref_allInfo=pd.DataFrame(columns=['Name','Wavelength', 'Intensity', 'Pixel'])
        for i in range(len(ref_sorted_peaks)):
            ref_allInfo=pd.concat([ref_allInfo,ref_sorted_peaks[i]],sort=False)
            
        return Material_names_sorted,ref_sorted_peaks,ref_allInfo

    ################################################## SECOND AVANTES SPECTROMETER RANGE
    
    elif Avantes==2:

        location=r'\\ipmlan\PK\Projekte\Fraunhofer\600610_InlineElement\10_akt_Arbeiten_Mitarbeiter\Alireza\Automatic libs sample analysis with database\505_740nm'
        X_train, _, _, Material_names_sorted,wavelengths=create_sorted_training_set(location,
                                                                                      peak_threshold=peak_threshold,
                                                                                      lam=lam,
                                                                                      positive=positive,
                                                                                      niteration=niteration)

        ref_sorted_peaks=[]
        for i in range(len(X_train)):
            idx_sorted,sorted_peaks=find_libs_sorted_peaks(X_train[i],
                                                           wavelengths,
                                                           peak_threshold=peak_threshold)
            ref_sorted_peaks.append(sorted_peaks)
            ref_sorted_peaks[i]['Name']=Material_names_sorted[i]

        ref_allInfo=pd.DataFrame(columns=['Name','Wavelength', 'Intensity', 'Pixel'])
        for i in range(len(ref_sorted_peaks)):
            ref_allInfo=pd.concat([ref_allInfo,ref_sorted_peaks[i]],sort=False)
            
        return Material_names_sorted,ref_sorted_peaks,ref_allInfo
    


# In[11]:


#Avantes1 270-540nm
Material_names_sorted1,ref_sorted_peaks1,ref_allInfo1=Analyse_And_Read_Updated_Database(
                                      Avantes=1,
                                      peak_threshold=150,
                                      lam=0.001,
                                      positive=0.000001,
                                      niteration=20)
#Avantes2 505_740nm
Material_names_sorted2,ref_sorted_peaks2,ref_allInfo2=Analyse_And_Read_Updated_Database(
                                      Avantes=2,
                                      peak_threshold=150,
                                      lam=0.001,
                                      positive=0.000001,
                                      niteration=20)


# In[12]:


def sample_refrence_comparison(spectrum,wavelengths,pixelshift=0,peak_threshold=150):

    idx_sorted,s_peaks_info=find_libs_sorted_peaks(spectrum,
                                        wavelengths,
                                        peak_threshold=peak_threshold)
    if wavelengths[0]==268.175:
        Material_names_sorted=Material_names_sorted1
        ref_sorted_peaks=ref_sorted_peaks1
        ref_allInfo=ref_allInfo1
    else:
        Material_names_sorted=Material_names_sorted2
        ref_sorted_peaks=ref_sorted_peaks2
        ref_allInfo=ref_allInfo2
    
    
    info=[]
    for s_peak_num in range(len(s_peaks_info)):

        s_peak_pixel=s_peaks_info.iloc[s_peak_num,2]
        s_peak_wavelength=s_peaks_info.iloc[s_peak_num,0]
        
        ###############cross refrencing
        ref_matching_info=pd.DataFrame(columns=['Name','Wavelength', 'Intensity', 'Pixel'])

        for material_num in range(len(ref_sorted_peaks)):

            ref_pixels=ref_sorted_peaks[material_num]['Pixel'].values

            ref_idx=np.where((ref_pixels==s_peak_pixel) 
                     | (ref_pixels==s_peak_pixel+pixelshift) 
                     | (ref_pixels==s_peak_pixel-pixelshift) )[0]

            if len(ref_idx) is not 0:

                ref_matching_info=pd.concat([ref_matching_info,ref_sorted_peaks[material_num].iloc[ref_idx]],
                                             sort=False)
                ref_matching_info=ref_matching_info.sort_values(by=['Intensity'],ascending=False)
                ref_matching_info.index

        info.append(ref_matching_info)
        info[s_peak_num].index=np.ones(len(info[s_peak_num]),'int8')*(int(s_peak_num)+1)

    #ALLinfo Dataframe
    allInfo=pd.DataFrame(columns=['Name','Wavelength', 'Intensity', 'Pixel'])
    for i in range(len(info)):
        allInfo=pd.concat([allInfo,info[i]],sort=False)

    #Final table Dataframe
    matched_mats=np.unique(allInfo['Name'].values)
    portions=[]
    s_portions=[]
    S_numOfPeaks=len(info)
    numOfmatchedMat=len(matched_mats)
    for i in range(len(matched_mats)):
        n_matched_peaks_for_mati=len(np.where(allInfo['Name'].values==matched_mats[i])[0])
        n_peaks_for_mati_in_ref=len(np.where(ref_allInfo['Name'].values==matched_mats[i])[0])
        portion=n_matched_peaks_for_mati/n_peaks_for_mati_in_ref
        s_portion=n_matched_peaks_for_mati/S_numOfPeaks
        portions.append(portion)
        s_portions.append(s_portion)

    portions_idX=np.asarray(portions).argsort()[::-1]
    sorted_portions=np.asarray(portions)[portions_idX]
    s_portions=np.asarray(s_portions)[portions_idX]
    sorted_mats=matched_mats[portions_idX]
    portions_info=[sorted_mat+ '| %0.2f'%sorted_portions[i]+':'+ '%0.2f'%s_portions[i]  for i,sorted_mat in enumerate(sorted_mats)]
    portions_info=pd.DataFrame(portions_info).T
    width=np.ceil(np.sqrt(numOfmatchedMat)).astype(int)
    final_table=np.chararray((width**2),unicode=True,itemsize=20)
    final_table[:portions_info.size]=portions_info
    final_table=pd.DataFrame(final_table.reshape(width,width))

    # Table Dataframe
    tableWidth=max([len(inf) for i , inf in enumerate(info)])
    table=np.chararray((len(info),tableWidth),unicode=True,itemsize=20)

    for i in range(len(info)):
        row=[]
        for i2 in range(len(info[i])):
            material=info[i].iloc[i2][0]
            ints=info[i].iloc[i2][2]
            row.append(material+':'+'%0.1f'%ints)
        table[i,:len(row)]=np.asarray(row)
    table=pd.DataFrame(table)
    table.index=np.arange(1,len(info)+1)

    return table, final_table,s_peaks_info,idx_sorted


# In[13]:


def full_info_plot(spectrum_wBG,spectrum_noBG,wavelengths,name,peak_threshold=150,pixelshift=0):

    table, final_table,s_peaks_info,idx_sorted= sample_refrence_comparison(spectrum_noBG,
                                                              wavelengths,
                                                              pixelshift=pixelshift,
                                                              peak_threshold=peak_threshold)
    #sat_idx=np.where(spectrum_Wbg[idx_sorted]>16000)[0]
    #sat_idx=np.asarray(sat_idx)
    
    ######## PLotting

    if len(s_peaks_info) is not 0:
            

        plt.figure(figsize=[100,50])   
        plotgrid=(15,21)

        ax=plt.subplot2grid(plotgrid, (1,3), colspan=5, rowspan=2)
        plt.plot(wavelengths,spectrum_wBG,linewidth=1,label='Spectrum')
        plt.xticks(np.arange(min(wavelengths), max(wavelengths)+1, 10.0))
        plt.title(name)
        plt.scatter(wavelengths[idx_sorted],spectrum_wBG[idx_sorted],linewidth=1,label='Peaks',marker="v",color='m')
        for i, pixel in enumerate(idx_sorted):
            plt.annotate(str(i+1), (wavelengths[pixel]-0.8, spectrum_wBG[pixel]),
                         size=9)
        plt.legend()

        ax=plt.subplot2grid(plotgrid, (3,3), colspan=5, rowspan=2)
        plt.plot(wavelengths,spectrum_noBG,linewidth=1,label='Spectrum without baseline')
        ticks=plt.xticks(np.arange(min(wavelengths), max(wavelengths)+1, 10.0))
        plt.legend()
        plt.scatter(wavelengths[idx_sorted],spectrum_noBG[idx_sorted],linewidth=1,label='Peaks',marker="v",color='m')
        for i, pixel in enumerate(idx_sorted):
            plt.annotate(str(i+1), (wavelengths[pixel]-0.8, spectrum_noBG[pixel]*1.02),
                         size=9)
         #############################################################################   
        tableLength=len(table)
        goodhight=50
        rowspan=np.ceil((tableLength/50)).astype(int)
        if tableLength>goodhight :
            tableLength=np.ceil((len(table)/2)).astype(int)

        ax=plt.subplot2grid(plotgrid, (0,8), colspan=6, rowspan=rowspan)
        ax.axis('off')
        the_table=plt.table(cellText=table.iloc[:tableLength,:].values,

                            rowLabels=table.iloc[:tableLength,:].index,
                            loc='upper left',
                            cellLoc='left')
        the_table.set_fontsize(20)
        the_table.scale(1, 1.5)

        if tableLength>goodhight :
            ax=plt.subplot2grid(plotgrid, (0,15), colspan=6, rowspan=rowspan)
            ax.axis('off')
            the_table=plt.table(cellText=table.iloc[tableLength:,:].values,

                                rowLabels=table.iloc[tableLength:,:].index,
                                loc='upper right',
                                cellLoc='left')
            the_table.set_fontsize(20)
            the_table.scale(1, 1.5)

        #plt.title('Matching peaks from Database(material:intensity in database)')

        ax=plt.subplot2grid(plotgrid, (6,4), colspan=3, rowspan=1)
        ax.axis('off')
        the_table=plt.table(cellText=final_table.values,      loc='center',                   cellLoc='left')
        plt.title('Material | observed fraction of peaks from database : fraction of the peaks in sample blonging to this material')
        the_table.set_fontsize(15)
        the_table.scale(1, 1.5)
        #########################################################

        ax=plt.subplot2grid(plotgrid, (1,0), colspan=3, rowspan=1)
        ax.axis('off')
        row=[]
        w=[]
        valPerCol=75
        if len(idx_sorted)<valPerCol:
            numOfcolumns=1
        elif len(idx_sorted)>valPerCol and len(idx_sorted)<valPerCol*2:
            numOfcolumns=2
        elif len(idx_sorted)>valPerCol*2 and len(idx_sorted)<valPerCol*3:
            numOfcolumns=3
        elif len(idx_sorted)>valPerCol*3 and len(idx_sorted)<valPerCol*4:
            numOfcolumns=4
        elif len(idx_sorted)>valPerCol*4 and len(idx_sorted)<valPerCol*5:
            numOfcolumns=5
        elif len(idx_sorted)>valPerCol*5 and len(idx_sorted)<valPerCol*6:
            numOfcolumns=6
        else:
            numOfcolumns=7
        c=0
        for i,lam in enumerate(wavelengths[idx_sorted]):
            row.append(str(i+1)+' | '+str(lam)+' : '+str(idx_sorted[i]))
            c+=1
            if c==numOfcolumns:
                w.append(row)
                c=0
                row=[]
        df=pd.DataFrame(w)
        the_table=plt.table(cellText=df.values,
                            loc='upper center',
                            cellLoc='center')
        the_table.set_fontsize(40)
        the_table.scale(.1866*numOfcolumns, 1.5)
        #plt.title('Peaks information from sample')

        plt.tight_layout()
    else:
                print(name+' has no peaks with the threshold '+str(peak_threshold))
    return final_table,s_peaks_info,table


# In[14]:


def Save_Plots_and_plotsDATABASEInfo(txtNames,
                                     wavelengths,
                                     spectra,
                                     saving_location,
                                     peak_threshold=150,
                                     pixelshift=0,
                                     lam=0.001,
                                     positive=0.000001,
                                     niteration=20,
                                     dpi=200):   
    ###saving location
    os.chdir(saving_location)
    
    spectra_no_BG=baseline_removal(spectra,lam=lam, positive=positive, niteration=niteration)

    spectra=np.atleast_2d(spectra)
    spectra_no_BG=np.atleast_2d(spectra_no_BG)
    
    #defining saving destinations
    saving_folder_name2='Peaks Information'+'  (peak threshold='+str(peak_threshold)+'pixelshift='+str(pixelshift)+')'
    saving_folder_name1='Plots'+'  (peak_threshold='+str(peak_threshold)+' pixelshift='+str(pixelshift)+')'

    for i in range(len(txtNames)):

        final_table,s_peaks_info,table=full_info_plot(spectra[i],
                                                       spectra_no_BG[i],
                                                       wavelengths,
                                                       name=txtNames[i],
                                                       peak_threshold=peak_threshold,
                                                       pixelshift=pixelshift)

        if len(s_peaks_info) is not 0:

                #saving the plots
            if not os.path.exists(saving_folder_name1):
                os.makedirs(saving_folder_name1)
            plt.savefig(saving_folder_name1+'/'+ txtNames[i]+'.jpg',dpi=dpi)
            plt.close()


                #Saving the peaks information
            if s_peaks_info is not None:
                if not os.path.exists(saving_folder_name2):
                    os.makedirs(saving_folder_name2)             
                s_peaks_info.to_csv(saving_folder_name2+'/'+ txtNames[i]+'_Sample_peaks_information.csv',header=True,index=False)
                final_table.to_csv(saving_folder_name2+'/'+ txtNames[i]+'_Material Peaks fractions.csv',header=True,index=False)
                table.to_csv(saving_folder_name2+'/'+ txtNames[i]+'_Database Comparison.csv',header=True,index=False)

                sorted_peaks=None
        #return idx_sorted,sorted_peaks  
    pass


# In[15]:


def makePlots_FullInfo_InSubfolders(location,
                                    peak_threshold=150,
                                    pixelshift=1,
                                    lam=0.001,
                                    positive=0.000001,
                                    niteration=20,
                                    dpi=200):                 
    
    subfolders=[]
    for root, dirs, files in os.walk(location):
        for name in files:
            if name.endswith((".txt")):
                subfolders.append(root)
                break

    for i in range(len(subfolders)):
        txtNames, wavelengths, spectra=AllTXTfilesSpectrumMatrix(subfolders[i])
        
        Save_Plots_and_plotsDATABASEInfo(txtNames,
                                         wavelengths,
                                         spectra,
                                         saving_location=subfolders[i],
                                         peak_threshold=peak_threshold,
                                         pixelshift=pixelshift,
                                         dpi=dpi)                             
            
    pass

