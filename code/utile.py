import numpy
from numpy import *
from sympy import Matrix
from math import sqrt
from test_primalité import millerRabin

#===========================================
def completer_matrice (M):
    A = mat(zeros(len(array(M)[0])))
    At = transpose(mat(zeros(len(M))))
    while(len(M)<len(array(M)[0])):
        M = concatenate((M, A), axis=0)
    while(len(M)>len(array(M)[0])):
        M = concatenate((M, At), axis=1)
    return M
        
def extraction_ligne (M):
    M = completer_matrice(M)
    A = array(M)
    M = Matrix(A)
    Mt = M.transpose()
    res = []
    i=0
    while(len(res)<=1):
        X = Mt.nullspace()[i]
        for j in range(len(X)):
            if (int(X[j])%2==1):
                res += [j]
        i+=1
    return res

def get_uv (l, M, n):
    list_x = extraction_ligne(M)
    u = 1
    v = 1
    for i in list_x:
        u = u*l[i]
        v = v*(l[i]**2 - n)
    return [u, get_racine(v)]
    


#===========================================
def eratostene(X) :
    tab = [x for x in range (X+1)]
    tab[1]=0

    for i in range (2, int(sqrt(X+1)+1)):
        if(tab[i]!=1):
            for j in range (2*i, X+1, i) :
                tab[j] = 1

    return [p for p in range (2, X+1) if tab[p]!=1]

        
def get_puissance_decomposition (l, B) :
    list_p = eratostene(B)
    res = [[0 for n in range (len(list_p))] for i in range (len(l))]
   
    for j in range (len(l)):
        for x in range(len(list_p)):
            while(l[j]%list_p[x] == 0):
                res[j][x]+=1
                l[j]=int(l[j]/list_p[x])
    return res

def cree_mat(liste_Q, B):
    tab = get_puissance_decomposition (liste_Q, B)
    for i in range (len(tab)):
        for j in range (len(tab[i])):
            tab[i][j] = tab[i][j]%2
    M = mat(tab)
    return M


#===========================================
def get_racine(x):
    r = int(sqrt(x))
    if(r**2>x):
        for i in range (len(str(r**2-x)), -1, -1):
            while(((r-10**i)**2)>x):
                r -= 10**i
        r = r-1        
    elif(r**2<x):
        for i in range (len(str(x-r**2)), -1, -1):
            while(((r+10**i)**2)<x):
                r += 10**i
    return r

#====================================================

#Retourne un tableau contenant tous les nombres B-friable entre 1 et X
def friable(B, X):
    tab=[i for i in range(1,X+1)] #On créé un tableau contenant les entiers de 1 a X
    for i in range (1, X):
        if(tab[i]!=1 and millerRabin(tab[i]) and tab[i]<=B):	#Si l'entier est premier, différent de 1 et inférieur à B, alors c'est un facteurs premiers qui nous intéresse
            for j in range(2*i+1, X, i+1):	#On fait des pats de taille i+1
                tab[j]=math.trunc(tab[j]/tab[i]) 	#On divise chaque multiple de tab[i] par tab[i]
            tab[i]=1
            
    return [p+1 for p in range (0, X-1) if tab[p]==1] 	#On prend les indices du tableau pour lequel la valeur est 1, il s'ait des entiers B-friables


def find_racine(T, p):
    b = False
    res = []
    for i in range(len(T)) :
        if ( T[i]%p == 0 ):
            res.append(i)
            if (b):
                break
            b = True
    return res

def find_quatre_racine(T, p):
    b = 1
    res = []
    for i in range(len(T)) :
        if ( T[i]%p == 0 ):
            res.append(i)
            if (b==4):
                break
            b = b+1
    return res


def find_racine_bis(r, p, n,k):
    res = []
    for a in r:
        A =(int)( a+ 1 + get_racine(n))
        for b in range (1, p):
            if(((2*A*b)%p) == 1 or ((2*A*b)%p) == (1-p)):
                res.append(int((A+(n-A**2)*b-get_racine(n))%p**(k+1))-1)
    return res



def friable_bis(B, X, n):
    r = get_racine(n)
    T = [int(((r+i)**2-n)) for i in range (1, X)]
    P = [p for p in range (3,B+1) if millerRabin(p)]
    x = 1
    go = True
    puissance = 1
    
    for p in P :
        racine = find_racine(T, p)
        for i in racine :
            for j in range(i, len(T), p):
                T[j] = int(T[j]/p)
        racine_bis=racine        
        while (go):
            racine_bis = find_racine_bis(racine_bis, p, n,puissance)
            puissance += 1
            go = False
            for i in racine_bis :
                for j in range(i, len(T), p**puissance):
                    T[j] = int(T[j]/p)
                    go = True
        go = True
        puissance=1

    p=2
        
    if( ((r+1)**2-n)%2 == 0 ):
        for i in range (0, len(T), 2):
            T[i] = int(T[i]/2)
    else :
        for i in range (1, len(T), 2):
            T[i] = int(T[i]/2)

    if( n%4 == 3 ):
        return [r+t+1 for t in range (len(T)) if (T[t] == 1)]
    else :
        if(n%8==5):
            racine=find_racine(T,2)
            for i in racine :
                for j in range(i, len(T), 4):
                    T[j] = int(T[j]/2) 
            return [r+t+1 for t in range (len(T)) if (T[t] == 1)]
        else:
            racine=find_quatre_racine(T, 2)
            for i in racine :
                for j in range(i, len(T), 4):
                   T[j] = int(T[j]/2) 
           
            if( n%8 == 1 ):
                racine = find_quatre_racine(T, 2)
                for i in racine :
                    for j in range(i, len(T), 8):
                        if(T[j] != 1):
                            T[j] = int(T[j]/2)
                puissance =4
                
                while (go):
                    racine = find_quatre_racine(T, 2)
                    go = False
                    for i in racine :
                        for j in range(i, len(T), 2**puissance):
                            if(T[j]!=1):
                                T[j] = int(T[j]/2)
                            go = True
                    puissance += 1
            go = True
            puissance = 2
                        
    return [r+t+1 for t in range (len(T)) if (T[t] == 1)]


#========================================================
def lancement():
    n = 39335476910299
    B = 179
    A = 3*10**5
    print("Les entiers de la formes Q(X)=X²-39335476910299 qui sont 179-friable avec X compris entre "+str(get_racine(n))+" et "+str(n)+" sont:")
    print(friable_bis(B, A, n))
