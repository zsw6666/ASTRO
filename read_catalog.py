#this script is use to read catalog and fit the data
import os
import test
import numpy
import math
from scipy import integrate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
os.chdir('/home/wu/download/Shiwu_astro/Flashlight_catalog/Images (5)')
Ilist=os.listdir(os.getcwd())
print Ilist
catalog=open(Ilist[0])
RA=[]
DEC=[]#RA and DEC is used to store the coordinate in the catalog
for i in catalog:
        numx=i.find(',')
        RA_center=float(i[:numx])
        i=i[numx+1:]
        numy=i.find(',')
        DEC_center=float(i[:numy])
        RA.append(RA_center)
        DEC.append(DEC_center)
N=float(len(RA))#N is used in integral constraint for w(theta)
uni_ra=numpy.random.uniform(low=min(RA),high=max(RA),size=(1,50*len(RA)))
uni_dec=numpy.random.uniform(low=min(DEC),high=max(DEC),size=(1,50*len(DEC)))
rr_ra=uni_ra[0]
rr_dec=uni_dec[0]
dr_ra=[];dr_dec=[]
dr_ra.extend(rr_ra);dr_ra.extend(RA)
dr_dec.extend(rr_dec);dr_dec.extend(DEC)#dr_ra and dr_dec are use to store the random and data's coordinate we combine them together at first to caculate the total number in each bin which include DD,RR,and DR,then we use this result subtract the DD and RR finally there is only DR remained
deta_d=0.001
[dd_num,dd_distance,mi,ma]=test.d_pairs(RA,DEC,deta_d)
[rr_num,rr_distance]=test.r_pairs(rr_ra,rr_dec,deta_d,mi,ma)
[dr_num_0,dr_distance]=test.r_pairs(dr_ra,dr_dec,deta_d,mi,ma)
DD_num=dd_num
dr_num_0=list(dr_num_0)
dr_num=[dr_num_0[i]-dd_num[i]-rr_num[i] for i in range(len(dr_num_0))]
dd_num=[i/(len(RA)*(len(RA)-1)*0.5) for i in dd_num]
rr_num=[i/(len(rr_ra)*(len(rr_ra)-1)*0.5) for i in rr_num]
s=(len(dr_ra)*(len(dr_ra)-1)*0.5)-(len(RA)*(len(RA)-1)*0.5)-(len(rr_ra)*(len(rr_ra)-1)*0.5);dr_num=[i/s for i in dr_num]#we normalize the dd_num,rr_num and dr_num
ww=[(dd_num[j]+rr_num[j]-2*dr_num[j])/rr_num[j] for j in range(len(dd_num)-1)]
ddistance=[(dd_distance[j]+rr_distance[j]+dr_distance[j])/float(3) for j in range(len(dd_distance)-1)]
w=[]
distance=[]
DD=[]
RR=[]
w_omiga=0
for i in range(len(ww)):
    if ww[i]>=0:
        w_omiga=w_omiga+dd_num[i]*ww[i]
        w.append(ww[i])
        RR.append(rr_num[i])
        DD.append(DD_num[i])
        distance.append(ddistance[i])
    else:#ww is used to store w(theta) value in each bin,due to the width of each bin is small and data sample is not very large,some values are minus,we abandon these value
        continue
ic1=[RR[i]*w[i] for i in range(len(w))]
IC=sum(ic1)/sum(RR)#IC is integral constraints
w=[i+IC+(1/N) for i in w]
errorbar=[math.sqrt((((1+w[i])/(1+w_omiga))**2)*(1/float(DD[i]))) for i in range(len(w))]
def func0(r,A):
    return A*(r**-0.8)
[A,pcov]=curve_fit(func0,distance,w)
A=A[0]
distance_fit=numpy.arange(min(distance),max(distance),(max(distance)-min(distance))/500)
w_fit=[func0(i,A) for i in distance_fit]
c=299792458;det_z=0.027;H0=7e4;Hr=3.68;W0=0.3;z0=2.254;sigma_8=0.8;gamma=1.8
print 'A='+str(A)
def func(z):
    return 1/((((1+z)**3)+(W0**-1)-1)**0.5)
inte=integrate.fixed_quad(func,0,z0)
x=(c/H0)*(W0**-0.5)*inte[0]
P=(W0**0.5)*(((1+z0)**3)+(1/W0)-1)**0.5
r0=((c*A*det_z)/(H0*Hr*(x**-0.8)*P))**(1.0/1.8)
print 'r0='+str(r0)

plt.figure(1)
plt.scatter(RA,DEC,marker='o',color='r',label='data')
plt.scatter(rr_ra,rr_dec,marker='x',color='b',label='random')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('distribution')
plt.legend(loc='upper right')
plt.show()
plt.figure(2)
plt.scatter(dd_distance,dd_num,marker='x',color='b',label='DD')
plt.scatter(rr_distance,rr_num,marker='o',color='r',label='RR')
plt.scatter(dr_distance,dr_num,marker='.',color='g',label='DR')
plt.xlabel('distance')
plt.ylabel('number')
plt.title('distance-number')
plt.legend(loc='upper right')
plt.show()
plt.figure(3)
plt.errorbar(distance,w,marker='v',color='b',yerr=errorbar)
plt.scatter(distance_fit,w_fit,marker='.',color='r')
plt.xlabel('distance')
plt.ylabel('w')
plt.title('angular correlation function')
plt.show()
