
#definition de function
def f(x):
    a=x[0]*x[0]+x[1]*x[1]
    return a


    
def dis(Points):
    sumx=0
    sumy=0
    for point in range(len(Points)):
        sumx=sumx+abs(Points[point][0]-Points[0][0])
        sumy=sumy+abs(Points[point][1]-Points[0][1])
    sum=sumx+sumy
    return sum


#choisir 3 point initiales
x1=[2,1]
x2=[2,2]
x3=[1,2]

#calcul les valeurs des points
v1=f(x1)
v2=f(x2)
v3=f(x3)


#order les point x1 va etre plus grand et x3 est plus petit
a=max(v1,v2,v3)

Points=[x1,x2,x3]
Valeurs=[v1,v2,v3]
print(dis(Points))
boucle=0
while boucle<15:
    boucle=boucle+1
    for point1 in range(len(Valeurs)):
        for point2 in range(len(Valeurs)):
            if Valeurs[point2]< Valeurs[point1]:
                item=Points[point2]
                itemv=Valeurs[point2]
                Points[point2]=Points[point1]
                Valeurs[point2]=Valeurs[point1]
                Points[point1]=item
                Valeurs[point1]=itemv
    #calcul centroid de xi sauf le plus grand
    print(Points)
    print(Valeurs)
    print("----------------------------------------------------------------")
    sumx=0
    sumy=0
    for i in range(len(Points)-1):
        sumx=sumx+Points[i+1][0]
        sumy=sumy+Points[i+1][0]    
    xcen=[sumx/(len(Points)-1.0),sumy/(len(Points)-1.0)]
    print(xcen,'here')
    #reflecion alpha=0.5
    alpha=1
    xrx=xcen[0]+alpha*(xcen[0]-Points[0][0])
    xry=xcen[1]+alpha*(xcen[1]-Points[0][1])
    xr=[xrx,xry]
    xrv=f(xr);
    #print(xr,f(xr))
    if xrv>=Valeurs[-1] and xrv < Valeurs[0]:
        Points[0]=xr
        print("reflection")
    #expansion gama=1.5   
    elif xrv<Valeurs[-1]:
       gama=1.5
       xex=xcen[0]+alpha*(xr[0]-xcen[0])
       xey=xcen[1]+alpha*(xr[1]-xcen[1])
       xe=[xex,xey]
       xev=f(xe);
       if xev<xrv:
           Points[0]=xe
           Valeurs[0]=xev
           print("expansion")
       else:
           Points[0]=xr
           Valeurs[0]=xrv
           print("expansion")
    #contraction ro=0.3
    elif xrv >= Valeurs[1]:
        ro=0.3
        xcx=xcen[0]-ro*(xr[0]-Points[0][0])
        xcy=xcen[1]-ro*(xr[1]-Points[0][1])
        xc=[xcx,xcy]
        xcv=f(xc)
        if xcv<Valeurs[0]:
            Points[0]=xc
            Valeurs[0]=xcv
            print("contraction")
    #shrink deta=0.5
    else:
        print("kun")
        for i in range(len(Points)-1):
            print('shrink')
            Points[i][0]=Points[-1][0]+0.5*(Points[i][0]-Points[-1][0])
            Points[i][1]=Points[-1][1]+0.5*(Points[i][1]-Points[-1][1])
    Valeurs=[f(Points[0]),f(Points[1]),f(Points[2])]
    



        

print(Points)
print(dis(Points))
print(boucle)
