import pylab as pl
from numpy import pi

# Percus-Yevick Terms  (see for instance D. Henderson Condensed Matter Physics 2009, Vol. 12, No 2, pp. 127-135)
def c0(eta):
    return -(1.+2.*eta)**2/(1.-eta)**4
def c1(eta):
    return 6.*eta*(1.+eta*0.5)**2/(1.-eta)**4
def c3(eta):
    return eta*0.5*c0(eta)
def cc(r,eta):
    if r>1:
        return 0
    else:
        return c0(eta)+c1(eta)*r+c3(eta)*r**3
# Spherical Fourier Transforms (using the liquid isotropicity)
def spherical_FT(f,k,r,dr):
    ft=pl.zeros(len(k))
    for i in range(len(k)):
        ft[i]=4.*pi*pl.sum(r*pl.sin(k[i]*r)*f*dr)/k[i]
    return ft

def inverse_spherical_FT(ff,k,r,dk):
    ift=pl.zeros(len(r))
    for i in range(len(r)):
        ift[i]=pl.sum(k*pl.sin(k*r[i])*ff*dk)/r[i]/(2*pi**2)
    return ift



def PercusYevickHS(phi,plot=True,filename="g_of_r.txt",x=pl.arange(0,5,0.05)):
    # number density
    rho=6./pi*phi
    # getting the direct correlation function c(r) from the analytic Percus-Yevick solution
    # vectorizing the function
    c=pl.vectorize(cc)
    # space discretization
    dr=0.005
    r=pl.arange(1,1024*2+1,1 )*dr
    # reciprocal space discretization (highest available frequency)
    dk=1/r[-1]
    k=pl.arange(1,1024*2+1,1 )*dk
    # direct correlation function c(r)
    c_direct=c(r,phi)
    # getting the Fourier transform
    ft_c_direct=spherical_FT(c_direct, k,r,dr)
    # using the Ornstein-Zernike equation, getting the structure factor
    ft_h=ft_c_direct/(1.-rho*ft_c_direct)
    # inverse Fourier transform
    h=inverse_spherical_FT(ft_h, k,r,dk)
    # print h
    # # radial distribution function
    gg=h+1
    # clean the r<1 region
    g=pl.zeros(len(gg))
    g[r>=1]=gg[r>=1]
    # save the cleaned version
    pl.savetxt(filename, zip(r,g))
    from scipy.interpolate import InterpolatedUnivariateSpline
    spl=InterpolatedUnivariateSpline(r, g)
    # plots
    if plot:
        pl.plot(r,abs((g-1)*r))
        # pl.plot(r,g-1)
        pl.ylabel("|(g(r)-1)r|")
        pl.xlabel("r")


    # return abs((spl(x)-1)*x)
    return spl(x)
# call the function
# PercusYevickHS(0.2)
# PercusYevickHS(0.45)
# PercusYevickHS(0.50)
# PercusYevickHS(0.49)
# bin,r,g=pl.loadtxt("exp.txt", unpack=True)





last=200

def fitting(ri,sigma,Phi):
    rr=ri/sigma
    return pl.log(pl.absolute(PercusYevickHS(phi=Phi, plot=False,x=rr[:last])-1)*rr)

from scipy.optimize import curve_fit


r,g,d,id2=pl.loadtxt("g.txt", unpack=True)
popt,pcov=curve_fit(fitting, r[20:last], pl.log(pl.absolute(g[20:last]-g[last])), p0=(12.,0.5))

print popt, pcov
# pl.plot(r[:last],g[:last])
PercusYevickHS(popt[-1],plot=True )
pl.plot(r[:last]/popt[0],abs((g[:last]-g[last]) *r[:last]/popt[0]))


# print len(spl)
# pl.plot(spl)

# pl.plot(r,( abs(g-g[len(g)/2])*r ),'o', ms=3., mfc='g', mec='none', mew=0.5, alpha=0.9)
# pl.xlim(0,6)
# pl.ylim(-1,3)
pl.yscale('log')
# PercusYevickHS(0.3)
pl.show()