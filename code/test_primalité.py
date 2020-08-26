import math
import random

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
 
# =============================================================================
def lpowmod(a, b, n):
    """exponentiation modulaire: calcule (a**b)%n"""
    r = 1
    while b>0:
        if b&1==0:
            b = b>>1
        else:
            r = (r*a)%n
            b = (b-1)>>1
        a = (a*a)%n
    return r
 
# =============================================================================
#  Test de primalité probabiliste de Miller-Rabin issue de http://python.jpvweb.com/python/mesrecettespython/doku.php?id=est_premier
def _millerRabin(a, n):
    """Ne pas appeler directement (fonction utilitaire). Appeler millerRabin(n, k=20)"""
    # trouver s et d pour transformer n-1 en (2**s)*d
    d = n - 1
    s = 0
    while d % 2 == 0:
        d >>= 1
        s += 1
 
    # calculer l'exponentiation modulaire (a**d)%n
    apow = lpowmod(a,d,n) # =(a**d)%n
 
    # si (a**d) % n ==1 => n est probablement 1er
    if apow == 1:
        return True
 
    for r in range(0,s):
        # si a**(d*(2**r)) % n == (n-1) => n est probablement 1er
        if lpowmod(a,d,n) == n-1:
            return True
        d *= 2
 
    return False
 
# ========================
def millerRabin(n, k=20):
    """Test de primalité probabiliste de Miller-Rabin"""
    global petitspremiers
 
    # éliminer le cas des petits nombres <=1024
    if n<=1024:
        if n in petitspremiers:
            return True
        else:
            return False
 
    # éliminer le cas des nombres pairs qui ne peuvent pas être 1ers!
    if n & 1 == 0:
        return False
 
    # recommencer le test k fois: seul les nb ayant réussi k fois seront True
    for repete in range(0, k):
        # trouver un nombre au hasard entre 1 et n-1 (bornes inclues)
        a = random.randint(1, n-1)
        # si le test echoue une seule fois => n est composé
        if not _millerRabin(a, n):
            return False
    # n a réussi les k tests => il est probablement 1er
    return True






#méthode des divisions succéssives
def div_suc (x):
    global petitspremiers
    if x<=1024:
        if x in petitspremiers:
            return True
        else:
            return False
    for i in range (2, math.trunc(math.sqrt(x))+1):
        if (x%i==0):
            return False
    return True
