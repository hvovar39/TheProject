import time
from math import *
from test_primalité import *
from sympy import gcd
from utile import *

petitspremiers = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,
    79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,
    179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,
    269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,
    367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,
    461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,
    571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,
    661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,
    773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,
    883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,
    1009,1013,1019,1021]


#Retourne la factorisation de n avec la méthode des divisions succéssives
def div_suc_fact (n):
    global petitspremiers
    if (millerRabin(n)):	#si n est premier, on renvoie n
        return [n]
    else:
        for i in petitspremiers:	#On commence par tester si de petits entiers divisent n (entre 2 et 1021)
            if(n%i==0):
                return div_suc_fact(int(n/i))+[i]
        for i in range (1021, math.trunc(sqrt(n))) :	#Si l'on n'a pas trouvé les facteurs dans les petits permiers, on vérifie les entiers de 1021 à sqrt(n)
            if(n%i==0):
                return div_suc_fact(int(n/i))+[i]
                

#====================================================

#factorise le nombre n en produit de nombre premier avec la méthode fermat
def fermat_bis(n):
    if(millerRabin(n)):
       return [n]
    for i in petitspremiers:
        if(n%i==0):
            return fermat_bis((int)(n/i))+[i]
    
    tab = fermat(n)
    r = tab[0]
    s = tab[1]
#On test si r-s et r+s sont bien des nombres premiers, sinon on relance fermat dessus
    if (millerRabin(r-s) and millerRabin(r+s)):
        return [r-s, r+s]
    elif (millerRabin(r-s)):
        return [r-s]+fermat_bis(r+s)
    elif (millerRabin(r+s)):
        return fermat_bis(r-s)+[r+s]
    else:
        return fermat_bis(r-s)+fermat_bis(r+s)

#donne une factorisation en 2 entiers(non forcément premier de n
def fermat(n):
    r=get_racine(n)+1	#On part de sqrt(n)+1
    while (get_racine((r**2)-n)**2 != ((r**2)-n) ):	#tant que r²-n n'est pas un carrée, on rajoute 1 a r
        r=r+1
    s=int(sqrt(r*r-n))	#s prend la valeur r²-n (qui est un entier de part le while précédent)
    return [r, s]

#====================================================
def cree_liste(n, B):
    r = get_racine(n)
    X = 2000
    A = int(B/log(B))
    liste_x = []
    liste_Q = []
    liste_friable = friable_bis(B, X, n)
    X = X*10
    for x in liste_friable:
        liste_Q.append(x**2-n)
        liste_x.append(x)
        A -= 1
    return [liste_x, liste_Q]

def crible_quadratique(n):
    if(millerRabin(n)):
       return [n]
    for i in petitpremiers:
        if(n%i==0):
            return fermat_bis((int)(n/i))+[i]
    
    B = int(exp((1/2)*sqrt(log(n)*log(log(n)))))
    listes = cree_liste(n, B)
    print(listes)
    M = cree_mat(listes[1], B)
    uv = get_uv(listes[0], M, n)
    u = uv[0]%n
    v = uv[1]%n
    p = gcd(max(u, v)-min(u, v), n)
    q = int(n/p)
    print(p)
    print(q)
    if(millerRabin(p) and millerRabin(q)):
        return [p, q]
    elif(millerRabin(p)):
        return [p]+crible_quadratique(q)
    elif(millerRabin(q)):
        return crible_quadratique(p)+[q]
    return crible_quadratique(p)+crible_quadratique(q)
    
    
#====================================================

#affiche le résultat de la factorisation de n par la méthode funct, ainsi que le temps d'éxecution
def factorisation(funct, n):
    time1 = time.clock()
    res = funct(n)
    time2 = time.clock() - time1
    s = ""
    for i in range(len(res)-1):
        s = s+str(res[i])+"*"
    s = s+str(res[len(res)-1])
    resultat = 1
    for j in res :
        resultat = resultat*j
    print ("La méthode à mis "+str(time2)+" pour trouver que :\n"+ str(n) + " = " + s)
    
#====================================================

tests_div = [2041, 34624234323236231, 2**83-1, 1099998619700431613, 2**47-1, 4333801, 39335476910299]
tests_fermat = [2041, 31885723060410621201917245580581940084008709974122013337, 1099998619700431613,  2**47-1, 4333801, 39335476910299]
tests = [2041, 34624234323236231, 2**83-1, 1099998619700431613, 4333801, 39335476910299, 31885723060410621201917245580581940084008709974122013337, 2**47-1]

print("====================================================================\nAvec la méthode des DIVISIONS SUCCESSIVES\n")
for n in tests_div :
    factorisation(div_suc_fact, n)

print()
print("\n===================================================================\nAvec la méthode de FERMAT\n") 
for n in tests_fermat :
    factorisation(fermat_bis, n)

