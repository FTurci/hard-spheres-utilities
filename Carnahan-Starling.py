#!/usr/local/bin/python
import sys
pi=3.14159265359
def cs(eta):
	return (1+eta+eta**2-eta**3)/(1-eta)**3*eta*6./pi

if __name__ == "__main__":
	nopt=len(sys.argv)
	if nopt<2:
		print "\n!! Provide an input packing fraction !!"
	if nopt==2:
	    P=cs(float(sys.argv[1]))
	    print P
	if nopt==3:
		import pylab as pl
		etas=pl.linspace(float(sys.argv[1]),float(sys.argv[2]),100)
		pl.plot(etas, cs(etas))
		pl.show()

