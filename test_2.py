def d_pairs(x,y,deta):
    '''
    This function is used to counts the galaxy pairs in each distance bin
    (x,y): coordinate of galaxies
    deta: width of distance bin
    '''
    import numpy
    
    
    #this part is used to calculate distance between each galaxy pairs
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
        coordinate.pop(0) 
    distance_database.sort()
    
    
    #in this part we define the distance bin
    distance_bin_width=deta
    k=int((distance_database[-1]-distance_database[0])/deta)
    distance_x=[distance_database[0]+i*distance_bin_width for i in range(1,k+2)]
    
    
    #this part is used to counts the number of pairs in each bin
    num=[]
    i=0;j=0;s=0;n=0
    while i<len(distance_database)-1 and j<len(distance_x):
        if distance_x[j]<distance_database[i]:
            num.append(n-s)
            s=n
            j+=1
            i-=1
        else:
            n+=1
        i+=1
    num.append(n-s)
    num.append(n-s)
    
    
    return num,distance_x,distance_database[0],distance_database[-1]

def r_pairs(x,y,deta,mi,ma):
    '''
    This function is used to statistic the random points pairs in each distance bin
    because we calculate two-point correlation function in the same distance bin,so,
    we need the same parameter who define bins with dd pairs
    (x,y): the coordinate
    deta,mi,ma: the parameter used to define distance bin
    '''
    import numpy



    #this part is used to calculate distance between each random point pairs
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
        coordinate.pop(0)
    distance_database.sort()
    
    
    #this part is used to define distance bin
    k=int((ma-mi)/deta)
    distance_x=[mi+i*deta for i in range(1,k+2)]
    
    
    #this part is used to counts the number of pairs in each bin
    num=[]
    i=0;j=0;s=0;n=0
    while i<len(distance_database)-1 and j<len(distance_x):
        if distance_x[j]<distance_database[i]:
            num.append(n-s)
            s=n
            j+=1
            i-=1
        else:
            n+=1
        i+=1
    return num,distance_x


def c_pairs(x,y,n):
    '''
    this function is calculate the galaxy pairs' counts within a certain distance 
    (x,y): the coordinate
    n; the number of bins we defined
    '''
    import numpy



    #this part is used to calculate distance between each galaxy pairs
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
        coordinate.pop(0) 
    distance_database.sort()
    
    
    #this part is used to define bin
    ma=max(distance_database)
    mi=min(distance_database)
    deta=(ma-mi)/n
    distance_a=numpy.arange(mi,ma+2*deta,deta)
    
    
    #this part is used to counts the number of galaxy pairs in each bins
    i=0;j=0;n=0
    num=[]
    while i<len(distance_database) and j<len(distance_a):
        if distance_a[j]<distance_database[i]:
            num.append(n)
        else:
            i+=1
        n+=1
    return num,distance_a,mi,ma

def cr_pairs(x,y,n,mi,ma):
    '''
    the same with r_pairs but used to counts the number of pairs within a certain distance 
    '''
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
        coordinate.pop(0)
    distance_database.sort()
    ini=distance_database[0]
    en=distance_database[-1]
    deta=(ma-mi)/n
    distance_a=numpy.arange(mi,ma+2*deta,deta)
    i=0;j=0;n=0
    num=[]
    while i<len(distance_database) and j<len(distance_a):
        if distance_a[j]<distance_database[i]:
            num.append(n)
        else:
            i+=1
        n+=1
    return num,distance,ini,en
