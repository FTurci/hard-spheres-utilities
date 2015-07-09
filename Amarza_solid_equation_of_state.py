import pylab as pl


def find_solid_p(rho):
    sigma=1.
    p=pl.arange(10,20,0.01 )

    v0=sigma**3/pl.sqrt(2)
    v=1./rho

    y=(p*sigma**3)**-1


    coeffs=pl.zeros(7)

    coeffs[0]=19328.09
    coeffs[1]=-2609.26
    coeffs[2]=141.6000
    coeffs[3]=11.56350
    coeffs[4]=-1.807846
    coeffs[5]=3
    coeffs[6]=-v+v0

    def find_real_positive(roots):
        for root in roots:
            if root.imag==0 and root.real>0:
                return root.real
    # print coeffs
    roots=1./pl.roots(coeffs)
    # f=3-1.807846*y+11.56350*y**2+141.6000*y**3-2609.260*y**4+19328.09*y**5-(v-v0)/y
    # # print pl.polyval(coeffs,1./p)

    # polynomial= pl.poly1d(coeffs)


    P=find_real_positive(roots)
    return P

def cs(rho):
    eta=rho/6.*pl.pi
    return (1+eta+eta**2-eta**3)/(1-eta)**3*eta*6./pl.pi

rhofirst=1.105
print "Starting density",rhofirst
print "Starting pressure",find_solid_p(1.105)

compressionfactor=0.98
rhos=[rhofirst*compressionfactor**n for n in range(1,12)]
print "densities",rhos
ps=[find_solid_p(rho) for rho in rhos]
pf=[cs(rho) for rho in rhos]
print "packing fractions", pl.array(rhos)/6*pl.pi
print "solid pressures",ps
print "fluid pressures",pf

print find_solid_p(.555014632)




# # pl.plot(p,f)
# pl.plot(p,pl.zeros(len(f)))
# pl.show()