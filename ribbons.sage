"""
Helper functions and such can go here.
"""
def ribbonBuddy(par):
    """
    Return the unique partition mu such that par/mu is a ribbon"""
    return Partition([par[i+1]-1 for i in range(len(par)-1)])

def ribbonVal(ribbon, seq):
    """
    Compute the determinant of a JT matrix indexed by a ribbon shape
    """
    a = lambda i: 0 if (i<0 or i>=len(seq)) else seq[i]
    lam = ribbon[0]
    ell = lam.length()
    mu = list(ribbon[1]) + (ell-ribbon[1].length())*[0]
    return matrix(ZZ, ell, ell, lambda i,j: a(lam[i]-mu[j]-i+j)).det()

def altRibbonVal(ribbon, seq):
    """
    Compute the determinant of a JT matrix indexed by a ribbon shape via symmetric functions
    """
    a = lambda i: 0 if (i<0 or i>=len(seq)) else seq[i]
    lam = ribbon[0]
    mu = ribbon[1]
    f = s(lam).skew_by(s(mu))
    return h(f).map_item(lambda b,c: (Partition([]), c*reduce(lambda x,y:x*y, [a(b[i]) for i in range(b.length())])))


def makeSequence(m,n,num):
    """
    Construct a sequence of the form [1, n, m, ...., m, n, 1]
    where m is repeated num times.
    """
    return [1, n] + [m]*num + [n,1]


def testAllRibbons(seq, bound):
    """
    Given a sequence ``seq``, compute the value of all ribbons of size at most ``bound``.
    Returns the first ribbon which gives a negative value, or ``None`` if no such ribbon is found.
    """
    for n in range(1,bound+1):
        for comp in Compositions(n):
            ribbon = comp.to_skew_partition()
            val = ribbonVal(ribbon, seq)
            if val < 0:
                print("OH NO!!", ribbon, val)
                return ribbon
    return None