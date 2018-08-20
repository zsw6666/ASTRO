import os
import pdb
import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u,constants as const
from astropy.cosmology import FlatLambdaCDM as FlatCDM

cosmo=FlatCDM(H0=68*u.km/u.s/u.Mpc,Tcmb0=2.725*u.K,Om0=0.3)


def IO(path,filename):
    '''
    This function is used to read data fro catalog
    path: the location of the catalog
    filename: name off the catalog
    '''

    os.chdir(path)
    data=np.hsplit(np.genfromtxt(filename,skip_header=1),8) 
    coordinate=np.array([data[1],data[2]]) 
    return coordinate


def auto_pairs(x,y,z0,n):
    '''
    This function is used to statistic the number of galaxy-galaxy pairs in each bins
    x: ra of each galaxy
    y: dec of each galaxy
    z0: redshift of quasar 
    n: number of bins
    '''
    
    x,y=np.array(x),np.array(y) 
    distanceset=[] 
    for i in range(len(x)):
        x0,y0=x[i],y[i]
        x1,y1=np.delete(x,i),np.delete(y,i) 
        deta_x,deta_y=x1-x0,y1-y0
        l=np.sqrt((deta_x**2)+(deta_y**2))
        distanceset=np.append(distanceset,l) #for each galaxy we calculate the distances between it and the rest
    counts,bins_a=Bin_Counts(distanceset,n,len(distanceset),True)#counts number of galaxy in each bin
    bins=bins_a*cosmo.comoving_distance(z0).value 
    return counts/len(x),bins,bins_a


def cross_pairs(x1,y1,x2,y2,z0):
    x1,y1,x2,y2=np.array(x1),np.array(y1),np.array(x2),np.array(y2)
    pointset1,pointset2=np.c_[x1,y1],np.c_[x2,y2] 
    distanceset=np.sqrt(np.sum((pointset2-pointset1)**2,1))  
    return distanceset 


def Bin_Counts(distanceset,binss,l,d_c):
   
    counts,bins=np.histogram(distanceset,binss,normed=d_c)
    if d_c:
        cd=np.cumsum(counts*np.diff(bins))*l  
        return cd,bins
    else:
        return counts,bins



def ACF(files,z0,n):
    '''
    This function is used to calculate angular correlation function
    and estimate mass of host halo by comparing bias we calculate and bias from simulation

    files:path+filesname
    z0:redshift of the target
    n: number of bins
    '''
   
    import halomod as hm
    from scipy.special import hyp2f1 as ghf
    from scipy.integrate import quad
    from scipy.optimize import curve_fit as cf 
    from scipy.optimize import leastsq as ltq
    
    #read data and create uniform field
    coordinate=IO(files[0],files[1])
    ra_g,dec_g=coordinate[0][1:]*np.pi/180.,coordinate[1][1:]*np.pi/180.
    ra_q,dec_q=coordinate[0][0]*np.pi/180.,coordinate[1][0]*np.pi/180.
    ra_r=np.random.uniform(low=np.min(ra_g),high=np.max(ra_g),size=(1,len(ra_g)))[0]
    dec_r=np.random.uniform(low=np.min(dec_g),high=np.max(dec_g),size=(1,len(dec_g)))[0] 
   

    #calculate cumulative galaxy-galaxy pairs and random-random pairs
    pairs_gg,bins_gg,bins_a=auto_pairs(ra_g,dec_g,z0,10)
    pairs_rr,bins_rr,bins_a=auto_pairs(ra_r,dec_r,z0,bins_a)
    A=np.pi*bins_gg[1:]**2
    ro_gg,ro_rr=pairs_gg/A,pairs_rr/A    


    #because there is statistical bias when counts pairs for different points for both galaxy-galaxy pairs and
    #random-random pairs,so,we divide ro_gg with ro_fit instead of ro_rr to get rid of this bias
    index_auto=np.where(bins_gg>=5)
    fit_func=np.polyfit(bins_rr[index_auto],ro_rr[(index_auto[0]-1)],1)
    ro_fit=fit_func[0]*bins_rr[1:]+fit_func[1]  
    index_auto=np.where(bins_gg<=4.5)
    bins_auto,ro_gg,ro_fit=bins_gg[(index_auto[0]+1)],ro_gg[index_auto],ro_fit[index_auto]
    w_auto=(ro_gg/ro_fit)-1 
    
    
    #calculate cumulative galaxy-qso pairs and qso-random pairs
    pairs_gq,bins_gq,bins_a=cross_pairs(ra_q,dec_q,ra_g,dec_g,z0,10) 
    A2=np.pi*bins_gq**2
    ro_gq=pairs_gq/A2
    ro_rq=len(ra_g)/A0
    index_cross=np.where(bins_gq<=2.2)
    bins_cross,ro_gq=bins_gq[index_cross],ro_gq[index_cross]
    w_cross=(ro_gq/ro_rq)-1 
   
   

    #define the fit function and fit the data
    def fitfunc(R,r0):
        gamma,deta_z=2.4,18
        a=((r0/R)**gamma)*ghf(0.5,0.5*gamma,1.5,-(deta_z/R)**2)
        return((r0/R)**gamma)*ghf(0.5,0.5*gamma,1.5,-(deta_z/R)**2)   
    def err(p,x,y):
        return fitfunc(x,p)-y
    corr_length_cross,pocov=cf(fitfunc,bins_cross,w_cross) 
    corr_length_auto=ltq(err,5,args=(bins_auto,w_auto))[0]    
  
    
    
    plt.figure('cross corr 1')
    plt.title('GQ density-radius relation')
    plt.xlabel('r Mpc')
    plt.ylabel('ro Mpc-2')
    plt.scatter(bins_cross,ro_gq,label='GQ')
    plt.scatter(bins_cross,np.full((1,len(bins_cross)),ro_rq),c='g',label='RQ')
    plt.legend()
    plt.figure('cross corr 2')
    plt.title('GQ cross correlation function')
    plt.xlabel('r Mpc')
    plt.ylabel('wp')
    plt.scatter(bins_cross,w_cross)
    plt.plot(bins_cross,fitfunc(bins_cross,*corr_length_cross),c='r')
    plt.figure('auto corr 1')
    plt.title('GG density-radius relation')
    plt.xlabel('r Mpc')
    plt.ylabel('ro Mpc-2')
    plt.scatter(bins_auto,ro_gg,label='GG')
    plt.scatter(bins_auto,ro_rr[index_auto],c='r',label='RR')
    plt.scatter(bins_auto,ro_fit,c='g',label='RR_FIT')
    plt.legend()
    plt.figure('auto corr 2')
    plt.title('GG auto correlation function')
    plt.xlabel('r Mpc')
    plt.ylabel('wp')
    plt.scatter(bins_auto,w_auto)
    plt.plot(bins_auto,fitfunc(bins_auto,corr_length_auto),c='r')
    plt.show()

#if __name__=='__main__':
#    halo_mass,bias,r0=ACF(['/home/wu/ASTRO/GCA/Flashlight_catalog/dense','SDSSJ0938+0905_LAE_cat.txt'],z0=2.254,n=20) 
