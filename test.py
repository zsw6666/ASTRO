# because there are two different kind of correlation functions(differential correlation function and cumulative correlation function) and to the counts the pairs in data catalog,random catalog and between them,I define four function(they are similiar) to do this
# (1) d_pairs(x,y,deta):counts pairs in data catalog for differential correlation function
# (2) r_pairs(x,y,deta,mi,ma):counts pairs in random catalog for differential correlation function
# (3) c_pairs(x,y,n):counts pairs in data catalog for cumulative correlation function
# (4) cr_pairs(x,y,n,mi,ma):counts pairs in random catalog for cumulative correlation function
#they all output two results,one is num which is a list to store counts in different distance bin the other is distance to store distance 
# the d_pairs and r_pairs is not so effient maybe we can first sort the distance and divide it in different bin and use the mean of the distance in the same bin as its distance 


def d_pairs(x,y,deta):#statistic pairs with different distance,x y are coordinate and deta is width of bin
    #we caculate distance of all possible pairs and put them in each bin 
    import numpy
    coordinate=[]
    for i in range(len(x)):
        o=[x[i],y[i]]
        coordinate.append(o)
    distance_database=[]#this list use to store all distance
    while len(coordinate)>1:#in this loop we firstly caclulate the distance between the first point with the rest and then pop the first out so that the original second become the first until there is one point in coordinate
        x0=coordinate[0][0]
        y0=coordinate[0][1]
        for i in range(1,len(coordinate)):
            x=coordinate[i][0]
            y=coordinate[i][1]
            l=(((x0-x)**2)+((y0-y)**2))**0.5
            distance_database.append(l)
        coordinate.pop(0)#statistic all distance 
    ma=max(distance_database)
    mi=min(distance_database)
    num=numpy.zeros(int((ma-mi)/deta)+1)
    distance=numpy.zeros(int((ma-mi)/deta)+1)#find the maximal and minimal distance and we can define the distance bin
    for i in distance_database:
        p=mi+deta
        s=mi
        k=0
        while p<=ma+deta:
            if i>=s and i<p:
                num[k]=num[k]+1
                distance[k]=distance[k]+i
                break
            else:
                p=p+deta#statistic number of distance in different range
                s=s+deta
                k=k+1
    for i in range(len(num)):
        if num[i]!=0:
            distance[i]=distance[i]/num[i]
        else:
            distance[i]=mi+(i+0.5)*deta
    return num,distance,mi,ma

def r_pairs(x,y,deta,mi,ma):#statistic pairs with different distance,x is x-axis,y is y-axis
    #r_pairs is the same with d_pairs
    import numpy
    coordinate=[]
    for i in range(len(x)):
        o=[x[i],y[i]]
        coordinate.append(o)
    distance_database=[]
    
    while len(coordinate)>1:
        x0=coordinate[0][0]
        y0=coordinate[0][1]
        for i in range(1,len(coordinate)):
            x=coordinate[i][0]
            y=coordinate[i][1]
            l=(((x0-x)**2)+((y0-y)**2))**0.5
            distance_database.append(l)
        coordinate.pop(0)#statistic all distance 
    num=numpy.zeros(int((ma-mi)/deta)+1)
    distance=numpy.zeros(int((ma-mi)/deta)+1)
    for i in distance_database:
        p=mi+deta
        s=mi
        k=0
        while p<=ma+deta:
            if i>=s and i<p:
                num[k]=num[k]+1
                distance[k]=distance[k]+i
                break
            else:
                p=p+deta#statistic number of distance in different range
                s=s+deta
                k=k+1
    for i in range(len(num)):
        if num[i]!=0:
            distance[i]=distance[i]/num[i]
        else:
            distance[i]=mi+(i+0.5)*deta
    return num,distance


def c_pairs(x,y,n):#x,y is coordinate and n is the number of bin which is the same as deta above 
    #c_pairs firstly caculate all distance and then counts points who are within the distance 
    #I sort the distance at first and then compare the distance in distance_database with the distance in distance _a 
    import numpy
    coordinate=[]
    for i in range(len(x)):
        o=[x[i],y[i]]
        coordinate.append(o)
    distance_database=[]
    while len(coordinate)>1:
        x0=coordinate[0][0]
        y0=coordinate[0][1]
        for i in range(1,len(coordinate)):
            x=coordinate[i][0]
            y=coordinate[i][1]
            l=(((x0-x)**2)+((y0-y)**2))**0.5
            distance_database.append(l)
        coordinate.pop(0)#statistic all distance 
    distance_database.sort()#sort all distance and then compare them with distance bin defined
    ma=max(distance_database)
    mi=min(distance_database)
    deta=(ma-mi)/n
    distance_a=numpy.arange(mi,ma+2*deta,deta)
    distance_iter=iter(distance_a)#use iter to generate iteration 
    num=[];l0=distance_iter.next();distance=[];i=0;k=0#i is the subscript of distance_database,l0 is the upper limit of each bin,k is used to counts the points within l0
    while i<=len(distance_database)-1:#in this loop,we first compare the elements in distance_database with the value in distance_a,if the former is smaller than latter then 
        if (i!=len(distance_database) and i!=len(distance_database)-1 and distance_database[i]<=l0 and distance_database[i+1]>l0) or i==len(distance_database)-1:
            num.append(k+1)
            distance.append(l0)
        if distance_database[i]<=l0:
            k=k+1
            i=i+1
        else:
            try:
                l0=distance_iter.next()
            except StopIteration:
                break
            else:
                continue
    return num,distance,mi,ma

def cr_pairs(x,y,n,mi,ma):
    #cr_pairs is the same with c_pairs
    import numpy
    coordinate=[]
    for i in range(len(x)):
        o=[x[i],y[i]]
        coordinate.append(o)
    distance_database=[]
    while len(coordinate)>1:
        x0=coordinate[0][0]
        y0=coordinate[0][1]
        for i in range(1,len(coordinate)):
            x=coordinate[i][0]
            y=coordinate[i][1]
            l=(((x0-x)**2)+((y0-y)**2))**0.5
            distance_database.append(l)
        coordinate.pop(0)#statistic all distance 
    distance_database.sort()
    ini=distance_database[0]
    en=distance_database[-1]
    deta=(ma-mi)/n
    distance_a=numpy.arange(mi,ma+2*deta,deta)
    distance_iter=iter(distance_a)
    num=[];l0=distance_iter.next();distance=[];i=0;k=0
    while i<=len(distance_database)-1:
        if (i!=len(distance_database) and i!=len(distance_database)-1 and distance_database[i]<=l0 and distance_database[i+1]>l0) or i==len(distance_database)-1:
            num.append(k+1)
            distance.append(l0)
        if distance_database[i]<=l0:
            k=k+1
            i=i+1
        else:
            try:
                l0=distance_iter.next()
            except StopIteration:
                break
            else:
                continue
    return num,distance,ini,en
if __name__=='__main__':
    import matplotlib.pyplot as plt
    x_ini=range(4)
    y_ini=range(4)
    x=[]
    y=[]
    for i in x_ini:
        for j in y_ini:
            x.append(i)
            y.append(j)
    [num,distance,mi,ma]=d_pairs(x,y,0.1)
    [num2,distance2,ma,mi]=c_pairs(x,y,500)
    for i in range(len(num)):
        if num[i]>0:
            print num[i]
            print distance[i]
    print 'num:'
    print num
    print 'distance:'
    print distance
    print 'num2='
    print num2
    print 'distance2='
    print distance2
    print ma,mi
    plt.figure(1)
    plt.scatter(x,y,marker='o',color='r')
    plt.show()
    #this is use to test the functions defined above
