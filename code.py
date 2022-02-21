
def isprime(n):
    if n < 2:
        return False
    k = 2
    while k*k<=n:
        if n%k == 0:
            return False
        k += 1
    return True


def logP(l):
    logprod = 0 
    for p in range(1, 10*l):
        if isprime(p):
            logprod += l*math.log(1+1.0/p)+math.log((1+((1-1.0/p)/(1+1.0/p))**l)/2)
    return logprod

import mpmath
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


#bound on the probability that the product is less than delta
def logfirstbound(l, delta):
    return logP(l)-l*math.log(1.0/delta)+0.0502*l

#bound on the probability that the tail is more than delta
def logsecondbound(r, k, sigma, delta):
    N0 = 72185376951205
    return logS(r,k,sigma)+2*k*math.log(1.0/delta)-math.log(k*(2-sigma)-1)-(k*(2-sigma)-1)*math.log(N0)


def finalbound(r, k, sigma, delta, l):
    print logfirstbound(l, delta)
    print logsecondbound(r, k, sigma, delta)


finalbound(10000, 40, 1.41, 0.12, 1100)






