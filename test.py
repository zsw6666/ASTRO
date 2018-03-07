def pairs(x,y):#statistic pairs with different distance,x is x-axis,y is y-axis
    import numpy
    coordinate=[]
    for i in range(len(x)):
        o=[x[i],y[i]]
        coordinate.append(o)
    distance_database=[]

    while len(coordinate)>1:
        for i in range(1,len(coordinate)):
            x0=coordinate[0][0]
            y0=coordinate[0][0]
            x=coordinate[i][0]
            y=coordinate[i][0]
            l=(((x0-x)**2)+((y0-y)**2))**0.5
            distance_database.append(l)
        coordinate.pop(0)#statistic all distance 
    ma=max(distance_database)
    mi=min(distance_database)
    deta=float((ma-mi)/20)
    num=numpy.zeros(20)
    for i in distance_database:
        p=mi
        k=0
        while p<(ma-deta):
            if i>=p and i<p+deta:
                num[k]=num[k]+1
            p=p+deta#statistic number of distance in different range
            k=k+1
    
    distance=[]
    for i in range(20):
        distance.append(mi+deta*i)

    return num,distance

if __name__=='__main__':
    #x=[1,2,3,4,5,6,7,8,9,10]
    #y=[11,12,13,14,15,16,17,18,19,20]
    #print pairs(x,y)
    import numpy
    import matplotlib.pyplot as plt
    x_g=10*numpy.random.randn(1,1000)
    y_g=10*numpy.random.randn(1,1000)
    x_u=10*numpy.random.rand(1,1000)
    y_u=10*numpy.random.rand(1,1000)#generate coodinate
    #print x_g[0],y_g[0],x_u[0],y_u[0]
    [g_pairs_num,g_distance]=pairs(x_g[0],y_g[0])
    [u_pairs_num,u_distance]=pairs(x_u[0],y_u[0])
    #print len(g_pair_num),len(g_distance)
    plt.figure(1)
    plt.title('normal distribution')
    plt.scatter(x_g,y_g,marker='.')
    plt.figure(2)
    plt.title('uniform dsitribution')
    plt.scatter(x_u,y_u,marker='.')
    plt.figure('normal distribution')
    plt.scatter(g_distance,g_pairs_num,marker='x')
    plt.xlabel('distance')
    plt.ylabel('number')
    plt.title('normal distribution')
    plt.show()
    plt.figure('uniform distribution')
    plt.scatter(u_distance,u_pairs_num)
    plt.xlabel('distance')
    plt.ylabel('number')
    plt.title('uniform dsitribution')
    plt.show()
