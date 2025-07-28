'''
Sagemath module to compute the Hilbert–Kunz function and multiplicity of an An singularity via the Han–Monsky algorithm.

I will keep extending the function to cover more diagonal hypersurfaces, although the algorithm is computationally heavy in many cases.

todo: add explanation of how to know when it will take time
'''
def lam(p,k):
    a = (p-1)/2
    if p%2 == 0 or k>p-1:
        return 'The first argument needs to be odd and the second must be strictly smaller.'
    l = []
    for i in range(p):
        for j in range(p):
            if j>=i and i+j<=p-1:
                i0 = i
                j0 = j
                x = lambda k: 1 if j0-i0<=k<=i0+j0 else 0
                lij = vector([x(k) for k in range(p)])
                l.append(lij)
            elif j-i >= 0:
                i0 = p-1-j
                j0 = p-1-i
                x = lambda k: 1 if j0-i0<=k<=i0+j0 else 0
                lij = vector([x(k) for k in range(p)])
                l.append(lij)
            else:
                j0 = i
                i0 = j
                if i0+j0<=p-1:
                    x = lambda k: 1 if j0-i0<=k<=i0+j0 else 0
                    lij = vector([x(k) for k in range(p)])
                    l.append(lij)
                else:
                    i1 = p-1-j0
                    j1 = p-1-i0
                    x = lambda k: 1 if j1-i1<=k<=i1+j1 else 0
                    lij = vector([x(k) for k in range(p)])
                    l.append(lij)
    l = [matrix([l[i+p*k] for i in range(p)]).transpose() for k in range(p)]
    return l[k]


def theta(p,M0):
    '''
    This function takes lambda_k in Lambda_n (as the matrix M0) and outputs the matrix of theta(lambda_k)
    in Lambda_n+1.
    '''
    n = int(log(M0.nrows())/log(p))
    error = 0
    k = list(M0[0]).index(1) #This is to detect which lambda the matrix M0 is representing. But this is easy.
    # It is just the location of the 1 in the first row/column. That's what this line does.
    # This is important to check parity on the next line:
    if k%2!=0: error = 1 # This is for parity.
    M = zero_matrix(p^(n+1),p^(n+1)) # I wouldn't say that these matrices are sparse...
    for a in range(p):
        for i in range(p^n): # We first give the values at elements lambda_ip
            for j in range(p^n):
                # lambda_kp lambda_ip = lambda_error_k * lambda_error_i * theta(lambda_k * lambda_i) =
                # = lambda_error_k * lambda_error_i * sum theta(lambda_j) 
                # = lambda_error_k * lambda_error_i * sum lambda_(pj+error_j)
                if i%2 == 0:
                    if j%2 == 0:
                        M[p*j+error*(p-1)+(-1)^error*a,p*i+a] = M0[j,i]
                    else:
                        M[p*j+p-1-error*(p-1)+(-1)^(1-error)*a,p*i+a] = M0[j,i]
                else:
                    # i is odd now...
                    if j%2 == 0:
                        M[p*j+p-1-error*(p-1)+(-1)^(1-error)*a,p*i+a] = M0[j,i]
                    else:
                        M[p*j+error*(p-1)+(-1)^error*a,p*i+a] = M0[j,i]
    return M

def scaler(p,M0):
    # This function takes the matrix of lambda_k in Lambda_n and gives you the matrix of lambda_k in Lambda_n+1
    return identity_matrix(p).tensor_product(M0)

def p_adic(p,n,k): # digits in base p as list [a_0,a_1,...] for a_0+a_1p+...
    k_p = []
    while k > 0:
        k_p.append(k%p)
        k = k//p
    if n>len(k_p):
        k_p = k_p + [0 for _ in range(n-len(k_p))]
    return k_p
# I need to check the error adjusments like p-1-error*(p-1)+(-1)^(1-error)*a. But it seems experimentally fine
#print((theta(p,theta(p,lam(p,k)))*scaler(3,scaler(3,lam(p,k)))).str(),'\n',scaler(3,scaler(3,lam(p,k))).str())

def lambda_matrix(p,n,k):
    if k>=p^n:
        return 'lambda_'+str(k)+'is not in Lambda_'+str(n)
    p_ad = p_adic(p,n,k)
    M = identity_matrix(p^n)
    for i in range(n): #len(p_ad) should be n!!
        Mi = lam(p,p_ad[i])
        for _ in range(i):
            Mi = theta(p,Mi)
        for _ in range(n-i-1):
            Mi = scaler(p,Mi)
        M = M*Mi
    return M

def vertical_swap(M):
    a = M.nrows()
    return matrix([[M[a-j-1,i] for i in range(a)] for j in range(a)])

def N_a(p,n):
    a = int((p^n-1)/2)
    return sum((-1)^k*2*lambda_matrix(p,n,k) for k in range(a))+(-1)^a*lambda_matrix(p,n,a)
def N_a_faster(p,n):
    '''
    I already know that this has a formula, so I should do a function that builds the matrix
    directly, and not from definition!
    '''
    return

def find_mu_e0(p,n,dim):
    '''
    The Han–Monsky's algorithm takes the exponents of a diagonal hypersurface, the characteristic,
    and a certain mu. If the correct mu is chosen, then the algorithm yields an expression for the
    Hilbert–Kunz function of the diagonal hypersurface that works for every e≥e_0 for a certain e_0.
    
    It can be seen that n+1 = c*p^e0 where gcd(c,p) = 1.
    
    This function computes the least mu such that p^mu = ± 1 (mod c). (see [HM93], Theorem 5.2)
    '''
    e0 = 0
    mu = 1
    while (n+1)%p^(e0+1) == 0: # This part of the function computes e0
        # e0 is the largest e such that p^e divides n+1
        # print((n+1)%p^(e0+1))
        e0 += 1
    c = (n+1)/p^e0
    if c!=1:
        #I think mu can be computed just thinking about the multiplicative group of units of Z/cZ
        while p^mu%c != 1 and p^mu%c != c-1%c: # But this computes it trying one by one
            # print(p^mu%c)
            mu+=1
    return mu,e0,c
    
def eHKAn(p,n,dim):
    '''
    Once you find the correct e0 and mu, then to compute the Hilbert–Kunz function one needs to
    use matrices and apply a formula.
    '''
    if n==1 and p == 2:
        m = dim//2
        if dim % 2 == 0:
            return (2**m+1)/2**m
        else:
            return 2**m/(2**m-1)
    elif p==2:
        raise ValueError('This function cannot compute examples in characteristic 2 (yet)')
    mu,e0,c = find_mu_e0(p,n,dim) 
    # n+1 = cp^e0
    # e0 is a shift needed in case n+1 is multiple of p
    # mu is the period
    beta = p^e0/(n+1)
    z = beta-floor(beta)
    z1 = p^mu*beta-floor(p^mu*beta)
    a = int((p^mu-1)/2)
    b = 0
    if c!=1:
        b = floor(p^(mu+e0)/(n+1))
    # print(mu,e0,c,z,z1,a,b)
    lsh = (lambda_matrix(p,mu,a)^dim*lambda_matrix(p,mu,b))[0,0]
    N_b = sum((-1)^k*c*lambda_matrix(p,mu+e0,k) for k in range(floor(p^mu/c)))+(-1)^(floor(p^mu/c))*(p^mu%c)*lambda_matrix(p,mu+e0,floor(p^mu/c))
    HK1 = p^e0*(N_a(p,e0)^dim)[0,0] # Ok
    HK2 = p^e0*((N_a(p,e0+mu)^dim)*N_b)[0,0] #If mu is greater than 3, this is already too big for p≥7.
    return (HK2 - HK1*lsh)/(p^(mu*dim)-lsh)/p^(e0*dim)#,(HK1)/(p^((e0)*dim))