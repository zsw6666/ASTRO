class BinaryTree:
    def __init__(self,rootobj):
        self.key=rootobj
        self.leftchild=None
        self.rightchild=None
    def insertleft(self,new_node):
        if self.leftchild==None:
            self.leftchild=new_node#the new_node given must be a tree 
                #we cannot use self.leftchild=new_node here because we must insert a object
        else:
            self.leftchild.insertleft(new_node)#if the left child is not empty,we keep searching the left child of this root until the left child is empty
    #insert right child like what we do above
    def insertright(self,new_node):
        if self.rightchild==None:
            self.rightchild=new_node
        else:
            self.rightchild.insertright(new_node)
    #visit root left child and right child
    def getleftchild(self):
        return self.leftchild

    def getrightchild(self):
        return self.rightchild

    def getroot(self):
        return self.key
    #set root's value
    def setrootvalue(self,value):
        self,key=value

def CreateKDTree(dot_database):#dot_database represent point in figure
    import numpy
    dot_database=numpy.array(dot_database)
    if len(dot_database)>1:
        xvar,yvar=numpy.var(dot_database,axis=0)#compare variance between x and y
        if yvar>=xvar:
            median=getmedian(dot_database,axis=1)#search for median
            root=getpoint(dot_database,median,axis=1)
            KDTree=BinaryTree(root)#build KDTree
            former,latter=get_subdatabase(dot_database,root,axis=1)#produce left-subdatabase and right-subdatabase
        else:
            median=getmedian(dot_database,axis=0)
            root=getpoint(dot_database,median,axis=0)
            KDTree=BinaryTree(root)
            former,latter=get_subdatabase(dot_database,root,axis=0)
        KDTree.insertleft(CreateKDTree(former))#insert leftchild and CreateKDTree produce a tree
        KDTree.insertright(CreateKDTree(latter))
        #create new KDTree instead of the original kdtree
    elif len(dot_database)==1:
        #print dot_database
        KDTree=BinaryTree(dot_database[0])
        root=dot_database
    else:
        KDTree=None
    return KDTree

def getmedian(database,axis):#obtain the median 
    if axis==1:
        y=[]
        for i in database:
            y.append(i[1])
        y.sort()
        median=y[int(len(y)/2)]
        return median
    else:
        x=[]
        for i in database:
            x.append(i[0])
        x.sort()
        median=x[int(len(x)/2)]
        return median


def getpoint(database,median,axis):#obtain the coodinate who own median caculated above
    if axis==1:
        for i in database:
            if i[1]==median:
                return i
    else:
        for i in database:
            if i[0]==median:
                return i

def get_subdatabase(database,root,axis):#produce the two subdatabase
    former=[]
    latter=[]
    database=[x for x in database if (x!=root).all()]#delete root.there is no pop method for numpy array
    if axis==1:
        for i in database:
            if i[1]<=root[1]:
                former.append(i)
            else:
                latter.append(i)
    else:
        for i in database:
            if i[0]<=root[0]:
                former.append(i)
            else:
                latter.append(i)
    return former,latter

if __name__=='__main__':
    a=[[1,2],[3,4],[5,6],[7,8],[9,10],[11,12],[13,14],[15,16],[17,18]]
    KDTree=CreateKDTree(a)
    print KDTree.getroot()
    print KDTree.getleftchild().getroot()
    print KDTree.getleftchild().getleftchild().getroot()
    print KDTree.getleftchild().getrightchild().getroot()
    print KDTree.getrightchild().getleftchild().getroot()
    print KDTree.getrightchild().getrightchild()
