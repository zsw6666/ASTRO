import os
import pdb
import numpy as np
from astropy import units as u,constants as const
from astropy.cosmology import FlatLambdaCDM as FlatCDM


def d_pairs(x,y):
    '''
    Thsis function is used to obtain 2-d two-point correlation function
    (x,y) is the coordinate of each galaxy
    '''
    
    
    #this part is used to calculate distance between every galaxy pairs in data
    x,y=np.array(x),np.array(y)
    pointset=np.c_[x,y]
    distanceset=np.array([])
    while len(pointset)>1:
        for i in range(1,len(pointset)):
            l=np.linalg.norm(pointset[0]-pointset[i])
            distanceset=np.append(distanceset,l)
        pointset=np.delete(pointset,0,0)
    
    #this part is used to statistic counts in each bin
    counts,bins=np.histogram(distanceset,bins='auto')
    return counts,bins


def r_pairs(x,y,bins):
    '''
    This function is the same with d_pairs but need one more parameter
    bins: we need this parameter because we must calculate dd pairs,rr_pairs and dr_pairs in the same bin
    '''
    x,y=np.array(x),np.array(y)
    pointset=np.c_[x,y]
    distanceset=np.array([])
    while len(pointset)>1:
        for i in range(1,len(pointset)):
            l=np.linalg.norm(pointset[0]-pointset[i])
            distanceset=np.append(distanceset,l)
        pointset=np.delete(pointset,0,0)
    
    counts,bins=np.histogram(distanceset,bins)
    
    return counts,bins


def LF(file_name,z,deta_z,sigma,deta_apha,deta_sigma,lamda_NB,deta_lamda_NB):
    '''
    This function is used to obtain the galaxy luminosity function
    file_name: the file's name(include path) which contains data 
    z: the redshift of the target
    deta_z: length of redshift interval
    sigma: DEC of the target
    deta_apha,deta_sigma: length of angular interval
    lamda_NB: effective wavelength of the narrowband filter
    deta_lamda_NB: FWHM of the emission line

    we use z deta_z sigma,deta_sigma,deta_apha to calculate the volume of the field which is used to calculate the number density
    lamda_NB and deta_lamda_NB are used to calculate the emission's flux

    '''


    #define some cosntant and convert some parameter to proper units
    cosmo=FlatCDM(H0=68*u.km/u.s/u.Mpc,Tcmb0=2.725*u.K,Om0=0.3)
    c=const.c;dl=cosmo.luminosity_distance(z);dl2=cosmo.luminosity_distance(z+deta_z);lamda_NB=lamda_NB*u.AA;deta_lamda_NB=deta_lamda_NB*u.AA;sigma=(sigma*u.deg).to(u.rad);deta_apha=(deta_apha*u.deg).to(u.rad);deta_sigma=(deta_sigma*u.deg).to(u.rad);deta_DL=dl2-dl
    deta_v_NB=(c/lamda_NB**2)*deta_lamda_NB
    
    
    #read the file which contains the magnitude data
    data=np.genfromtxt(file_name,delimiter=',')
    m_NB=data[:,0];m_g=data[:,1]
    

    #calculate the luminosity of each galaxy,the equation you can found here:
    #https://github.com/zsw6666/test/blob/master/Equation_F_LF.jpg
    F_Lya=(3631*u.Jy)*deta_v_NB*((10**(m_NB/-2.5))-(10**(m_g/-2.5)))
    L_Lya=(F_Lya*4*np.pi*dl**2).to(u.erg/u.s)
    
    
    #statistic counts in each bin and convert the counts to number density
    [counts,bins]=np.histogram(L_Lya,bins='auto') 
    Vmax=(dl**2)*np.cos(sigma)*deta_apha*deta_sigma*deta_DL
    deta_L=bins[1]-bins[0]
    counts=counts/Vmax 
    bins=np.delete(bins,np.append(np.where(counts==0.),0))*(u.erg/u.s)
    counts=np.delete(counts,np.where(counts==0.))

    
    errorbar=np.sqrt(counts)*(1/np.sqrt(Vmax))
    
    return bins.value,counts.value,errorbar.value,deta_L

def SCF(L,Le,apha,phie):
    '''
    Schechter luminosity function which is used to fit data
    Le,apha,phie: the parameter(phie is the normalization)
    L: the variable
    '''

    phi=(1/Le)*phie*((L/Le)**apha)*np.exp(-L/Le)
    return phi


#Example
if __name__=='__main__':
    import pdb
    import matplotlib.pyplot as plt
    L_lya,n,err,deta_L=LF('/home/wu/ASTRO/GCA/data/Flashlight_catalog/Images (5)/LAE-Mab5.txt',2.255,0.027,0.90915,0.009155,0.0095589,3955,32.7)
    index=[len(L_lya)-1,len(L_lya)-2] 
    L_lya=np.delete(L_lya,index)
    n=np.delete(n,index)
    err=np.delete(err,index) 
    L_fit=np.arange(min(L_lya),max(L_lya)+0.5*max(L_lya),(max(L_lya)-min(L_lya))/500.0)
    Le=10**42.33;apha=-1.65;phie=10**(-2.86)
    N_fit=SCF(L_fit,Le,apha,phie)*deta_L
    #pdb.set_trace()
    plt.figure(1)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.title('Luminosity Function')
    plt.xlabel('Luminosity erg/s')
    plt.ylabel('n Mpc-3')
    plt.xlim(min(L_lya),max(L_lya))
    plt.ylim(min(N_fit),max(n))
    plt.errorbar(L_lya,n,fmt='o',yerr=err)
    plt.plot(L_fit,N_fit,c='r')
    plt.show()
