def LF(Path,file_name,dL1,dL2,sigma,deta_alpha,deta_sigma,z):
    #Path:the Path of file which contain the target's magnitude of broad band and narrow band
    #file_name:the file's name
    #dL1:luminosity distance 1 /Mpc
    #dL2:luminosity distance 2 /Mpc
    #sigma:DEC
    #deta_apha:deta RA
    #deta_sigma:deta DEC
    #z:redshift
    
    
    import os
    import numpy
    
    #we read the file containing magnitude
    os.chdir(Path)
    File=open(file_name)
    m_NB=[];m_g=[]
    for i in File:
        index=i.find(',')
        m_NB.append(float(i[:index]))
        m_g.append(float(i[index+1:-1]))
    m_NB.sort();m_g.sort()#sorting the magnitude for convenience
    m_NB=m_NB[1:];m_g=m_g[1:]
    
    #kk set the number of bins
    #after setting kk,we can obtain the width of the bin(bin_NB and bin_g) according to the following equation
    #M_g and M_NB store the value of each bin(magnitude)
    kk=12
    bin_NB=(m_NB[-1]-m_NB[0])/kk;bin_g=(m_g[-1]-m_g[0])/kk
    M_NB=[m_NB[0]+k*bin_NB for k in range(1,kk+1)]
    M_g=[m_g[0]+k*bin_g for k in range(1,kk+1)]
    
    #N_g and N_NB contain the counts of each bin
    #because we have sort the magnitude(M_g and M_NB),we will be able to counts the number of LAEs in each magnitude bin 
    #traversing M_g and M_NB once.
    #The time comlexity is o(nlogn) including the time we spend on sorting 
    N_NB=[];N_g=[]
    M=[M_NB,M_g]
    N=[N_NB,N_g]
    m=[m_NB,m_g]
    for j in range(len(M)):
        k=0;i=0;n=0;s=0
        while i<len(m_g)-1 and k<kk:
            if M[j][k]<m[j][i]:
                N[j].append(n-s)
                s=n
                k+=1
                i-=1
            else:
                n+=1
            i+=1
        N[j].append(n-s)
    Mg=M[1];MNB=M[0]


    #The following steps are converting the magnitude to luminosity
    #lamada1 and lamda2 is the limit of g band
    #deta_v is the corresponding frequency width of this band
    #F_g and L_g contain flux and luminosity of each LAE in broad band respectively
    #lamda_NB is the effective wavelength of the narrow band filter 
    #deta_lamda_NB is the FWHM of the emission line
    #F_NB and L_NB contain flux and luminosity of each LAE in narrow band respectively
    #f_continuum is the continuum component of galaxy emission
    
    #because there are two kind of component in the narrow band(emission line and continuum) we need to substract the continuum component
    #so,we firstly calculate the average flux density of g band and regard it as the continuum component.
    #then we convert the NB magnitude to flux and substract it with the the continuum component
    #we can then calculate the narrow band luminosity of LAEs with the luminosity distance
    #because we should know the galaxy number density,we then calculate the volume(Vmax) using deta_apha,deta_sigma,sigma,and deta_DL
    lamda1=3620e-10;lamda2=5620e-10;c=299792458;deta_v=(c/lamda1)-(c/lamda2)
    F_g=[(1e4)*(1e-23)*3631*deta_v*(10**(i/(-2.5))) for i in Mg]#unit erg/s/m^2.
    L_g=[(4*numpy.pi)*((dL1*3.08567758e22)**2)*i for i in F_g]#unit erg/s
    lamda_NB=(3955e-10);deta_lamda_NB=32.7e-10
    deta_v_NB=(c/lamda_NB**2)*deta_lamda_NB
    F_NB=[(1e4)*(1e-23)*deta_v_NB*3631*(10**(i/(-2.5))) for i in MNB]#unit erg/s/m^2
    f_continuum=[(F_g[i]-F_NB[i])/(deta_v) for i in range(len(N_g))]
    F_NB=[F_NB[i]-f_continuum[i]*deta_v_NB for i in range(len(F_NB))]
    L_NB=[(4*numpy.pi)*((dL1*3.08567758e22)**2)*i for i in F_NB]#unit erg/s
    DL=(dL1+dL2)/2;deta_DL=dL2-dL1
    Vmax=(DL**2)*numpy.cos(sigma)*deta_DL*deta_alpha*deta_sigma
    
    #because the LAE sample of this target is small(125 LAEs) there are not sufficient LAEs in some bin,we then delete this bin
    i=0
    while i<len(N_NB):
        if N_NB[i]<=1.0:
            N_NB.pop(i)
            L_NB.pop(i)
        else:
            i+=1
    errorbar_NB=[numpy.sqrt(i)/Vmax for i in N_NB]
    N_NB=[i/Vmax for i in N[0]]
    N_g=[i/Vmax for i in N[1]]
    return L_NB,L_g,N_NB,N_g,errorbar_NB
def SCF(L,Le,apha,phie):
    phi=phie*((L/Le)**apha)*numpy.exp(-L/Le)
    return phi
if __name__=='__main__':
    import numpy
    import matplotlib.pyplot as plt
    a=LF('/home/wu/ASTRO/GCA/data/Flashlight_catalog/Images (5)','LAE-Mab5.csv',18424.07,18689.29,0.90915*(numpy.pi/180.0),0.009155*(numpy.pi/180.0),0.0095589*(numpy.pi/180.0),2.254)
    L_NB=a[0]
    N_NB=a[2]
    errorbar_NB=a[4]

    L_fit=numpy.arange(min(L_NB),max(L_NB)+0.5*max(L_NB),(max(L_NB)-min(L_NB))/500.0)
    Le=10**42.33;apha=-1.65;phie=10**(-2.86)
    N_fit=[SCF(L,Le,apha,phie) for L in L_fit]


    plt.figure(1)
    ax=plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    plt.xlim(min(L_NB),max(L_NB))
    plt.ylim(min(N_fit),max(N_fit))
    plt.errorbar(L_NB,N_NB,fmt='o',yerr=errorbar_NB)
    plt.plot(L_fit,N_fit,c='r')
    plt.show()
