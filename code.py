import math
import mpmath

def isprime(n):
    if n < 2:
        return False
    k = 2
    while k*k<=n:
        if n%k == 0:
            return False
        k += 1
    return True

#computes log(P_lambda). Here l = lambda since lambda is not permitted as a variable name in Python.
def logP(l):
    logprod = 0 
    for p in range(1, 10*l):
        if isprime(p):
            logprod += l*math.log(1+1.0/p)+math.log((1+((1-1.0/p)/(1+1.0/p))**l)/2)
            #this twisted expresion of the terms in the product is meant to avoid dealing with the large number (1+1/p)^l before log.
    return logprod

def zeta(s):
    return float(mpmath.zeta(s))

#computes log(S(R,k,sigma))
def logS(r, k, sigma):
    logprod = 0
    for p in range(1, r+1):
        if isprime(p):
            logprod += math.log(((1-1.0/p**(sigma/2))**(-2*k)+(1+1.0/p**(sigma/2))**(-2*k))/2)
            logprod += (2*k*k+2*k)*math.log(1-1.0/p**sigma)
    logprod += (2*k*k+2*k)*math.log(zeta(sigma))
    return logprod


#bound on the probability that the euler product is less than delta
def logfirstbound(l, delta):
    return logP(l)-l*math.log(1.0/delta)+0.0502*l

#bound on the probability that the tail is more than delta
def logsecondbound(r, k, sigma, delta):
    N0 = 72185376951205
    return logS(r,k,sigma)+2*k*math.log(1.0/delta)-math.log(k*(2-sigma)-1)-(k*(2-sigma)-1)*math.log(N0)

#returns the final bound
#the parameters are optimized by balancing the two bounds - we did it by tweaking the parameters by hand until the two bounds
#were of similar magnitude
#we include the third term math.exp(-7000), even though python handles this negligible contribution as 0.0
def finalbound(r, k, sigma, delta, l):
    l1 = logfirstbound(l, delta)
    l2 = logsecondbound(r, k, sigma, delta)
    print l1, l2
    return math.exp(l1)+math.exp(l2)+3*math.exp(-7000)


#this performs a random descent starting at a givn choice - together with analysis by hand it is helpful in the optimization
def optimize(r, k, sigma, delta, l, lim):
    currentrecord = 100
    oldtuple = (k, sigma, delta, l)
    for i in range(lim):
        newtuple = (oldtuple[0]+random.randint(-1,1), oldtuple[1]+0.0001*(random.random()-0.5), oldtuple[2]+0.0001*(random.random()-0.5), oldtuple[3]+random.randint(-1,1))
        if finalbound(r, newtuple[0], newtuple[1], newtuple[2], newtuple[3]) < currentrecord:
            currentrecord = finalbound(r, newtuple[0], newtuple[1], newtuple[2], newtuple[3])
            oldtuple = newtuple
            print newtuple
            print currentrecord
    
#this is roughly optimal choice we found, via a random descent

finalbound(10000, 48, 1.42, 0.12, 700)














