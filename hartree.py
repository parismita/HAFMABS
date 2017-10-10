import copy
import sys
import math
import copy
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------------------------------------------------------#
'''
This software is licensed under the Apache 2 license, quoted below.

Copyright 2017, Parismita Das <das.parismita@gmail.com>
Copyright 2017, Hussain Rasiwala <hussainsr97@gmail.com>

Licensed under the Apache License, Version 2.0 (the "License"); you may not
use this file except in compliance with the License. You may obtain a copy of
the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
License for the specific language governing permissions and limitations under
the License.'''
#------------------------------------------------------------------------------------------------------------------------------#
'''
###ALL THE VALUES ARE IN Ry atomic UNIT###

Variables:
mesh: grid points number
r: radial distance r
dx: diff b/w 2 points
zmesh: atomic number zmesh
vhx: potential due to charge density
vpot,ppp: total potencial
ei,e: energy eigenvalue
y: wavefunction
rho: charge density function

Functions:
do_mesh : print grid information such as r, dx, r(0), r(end)
init_pot: initialize potential as -2z/r (in Ry unit)
rho_of_r: charge density calculation for every r
v_of_rho: potencial function integrating charge density to get the potencial at r
solve_sheq : solving shrodingers equation using numerov's method in radial coordinates on a logarithmic grid.

Main Code:
calculating zeroth order energy eigenvalue , writing wavefunction and energy to file.
reading the wavefunction from file itteratively
using it to calculate the potential function
putting the new potential function to function solve_sheq
getting next wavefunctions and energy from schroedinger's eqn

Note:   In the code to verify the method the ionisation energy values are calculated 
		assuming we are removing the electrons to get the 1st,2nd and 3rd IP

		ground state energies calculated considering effect of all elecrons 


Final energy values in (eV) unit.
'''
#initialises distance r, computes r_square, r_sqrt
def do_mesh(mesh, zmesh, xmin, dx, rmax, r, sqr, r2):
	#initialize grid
	for i in range(mesh+1):
		x = xmin + dx * i;
		r[i] = math.exp(x) / zmesh;
		sqr[i] = math.sqrt(r[i]);
		r2[i] = r[i] * r[i];

	print " radial grid information:"
	print " dx   = ", dx
	print " xmin = ", xmin
	print " zmesh =", zmesh
	print " mesh = ", mesh
	print " r(0) = ",  r[0]
	print " r(mesh) = ", r[mesh]

############################################################################################################
#initialise zeroth potential as -2z/r
def init_pot(zeta,mesh,r,vpot):
	#initialize potential
	out = open("pot.out","w")
	out.write("r                    V(r)")
	out.write("\n")
	for i in range(mesh+1):
		vpot[i] = -2 * (zeta / r[i])
		out.write(str(r[i])+" " +str(vpot[i]))
		#out.write(str(vpot[i]) + " ")
		out.write("\n")

###########################################################################################################
#calculating the charge density as func of r
def rho_of_r(mesh, r, r2, y, rho):					    #y needs to be read from a file	
	nelec=2.0			     			             #here i changed it from 2 to 1
	fpi=4.0*3.14159265358945
	for i in range(mesh+1):
		rho[i] = nelec *y[i]*y[i]* r[i] / (fpi*r2[i]);			#might need to remove nelec


############################################################################################################
#potential func calculated
def v_of_rho(mesh,dx,r,r2,y):
	rho=[0 for j in range(mesh+1)]
	vhx=[0 for j in range(mesh+1)]

	rho_of_r(mesh, r, r2, y, rho)		
	e2=1.0 						               
	charge = 0.0
	fpi = 4.0*3.14159265358945
	for i in range(mesh+1):					   
		charge = charge + rho[i] * fpi * r2[i] * r[i] * dx		      
    		vhx[i] = e2*charge/r2[i]				     # electric field from guass law	
	
	ehx = 0.0;
        for i in range(mesh-1,-1,-1):
        	vhx[i] = vhx[i+1] + vhx[i] * r[i] * dx;      #calculating the potential coming from infinity
  
        return vhx
############################################################################################################
def solve_sheq(n, l, zeta, mesh, dx, r, sqr, r2, vpot, y):
	eps=1e-10
	n_iter=100
	ddx12 = dx * dx / 12.
	sqlhf = (l + 0.5) * (l + 0.5)
	x2l2 = (2*l+ 2)
	
	#set initial lower and upper bounds to the eigenvalue
	eup = vpot[mesh]
	elw = eup
	for i in range(mesh+1):
		elw = min(elw, sqlhf/r2[i]+vpot[i])
	if (eup - elw < eps):
		return 0
	e = (elw + eup) * .5
	f = [0.0 for i in range(mesh+1)]
				
	#start loop on energy	
	de= 1e+10
	for kkk in range(n_iter):
		if(abs(de)<eps):
			break
		icl = -1;
		f[0] = ddx12 * (sqlhf + r2[0] * (vpot[0] - e))
		for i in range(mesh+1):
			f[i] = ddx12 * (sqlhf + r2[i] * (vpot[i] - e))
			if (f[i] == 0.):
				f[i] = 1e-20;
			if (f[i] != math.copysign(f[i], f[i - 1])):
				icl = i;
		if ((icl < 0 )or (icl >= mesh - 2)):
			print "error in Schrodinger: last change of sign too far"
			return 0

		#f function as required by numerov method
		for i in range(mesh+1):
			f[i] = 1. - f[i]
			y[i] = 0.0

		nodes = n - l - 1;
		y[0] = pow (r[0], l+1) * (1. - zeta * 2. * r[0] / x2l2) / sqr[0]
		y[1] = pow (r[1], l+1) * (1. - zeta * 2. * r[1] / x2l2) / sqr[1]

		#outward integration, count number of crossings
		ncross = 0;
		for i in range(1,icl):
			y[i + 1] = ((12. - f[i] * 10.) * y[i] - f[i - 1] * y[i - 1])/f[i + 1]
			if (y[i] != math.copysign(y[i],y[i+1]) ):
				ncross = ncross+1

		fac = y[icl]
		#check number of crossings
		#/incorrect number of nodes: adjust energy bounds
		if (ncross != nodes):
			if (ncross > nodes):
				eup = e;
			else :
				elw = e;
			e = (eup + elw) * .5;
		else:
			#correct number of nodes: perform inward iteration
			y[mesh] = dx
			y[mesh - 1] = (12. - f[mesh] * 10.) * y[mesh] / f[mesh - 1]

			#inward integration
			for i in range(mesh-1,icl,-1):
				y[i - 1] = ((12. - f[i] * 10.) * y[i] - f[i + 1] * y[i + 1])/ f[i - 1];
				if (y[i - 1] > 1e10):
					for j in range(mesh-1,i-2,-1):
						y[j] /= y[i - 1];

			#rescale function to match at the classical turning point (icl)
			fac /= y[icl]
			for i in range(icl,mesh+1):
				y[i] *= fac

			#normalize on the segment
			norm = 0.
			for i in range(1,mesh+1):
				norm += y[i] * y[i] * r2[i] * dx
			norm = math.sqrt(norm)
			for i in range(mesh+1):
				y[i] /= norm

			#find the value of the cusp at the matching point (icl)
			i = icl
			ycusp = (y[i - 1] * f[i - 1] + f[i + 1] * y[i + 1] + f[i] * 10. * y[i]) / 12.
			dfcusp = f[i] * (y[i] / ycusp - 1.)

			#eigenvalue update using perturbation theory
			de = dfcusp / ddx12 * ycusp * ycusp * dx;
			if (de > 0.) :
				elw = e
			if (de < 0.) :
				eup = e

			e = e + de
			e = min(e,eup)
			e = max(e,elw)

	#convergence not achieved 
	if ( abs(de) > eps ):
		print " Schrodinger not converged after", n_iter,"iterations\n"
		return 0

#	print "convergence achieved at iter # ",kkk, "de = ",de
	
	return e



############################################################################################################
#initialize atomic charge (Z)
#print "\nPlease enter atomic charge :"
#zeta = int(raw_input().strip())

#initialize logarithmic mesh
zeta=3
zmesh = zeta;
rmax = 100.
xmin = -8.
dx = 0.01

#number of grid points
mesh = int((math.log(zmesh * rmax) - xmin) / dx)
r = [0.0 for i in range(mesh+1)]
r2= [0.0 for i in range(mesh+1)]
sqr= [0.0 for i in range(mesh+1)]
n=[1,1,2]
l=0.0
y=[]

do_mesh(mesh, zmesh, xmin, dx, rmax, r, sqr, r2)

#initialize the potential
vpot= [0.0 for i in range(mesh+1)]
init_pot(zeta, mesh, r, vpot);

#read number of nodes (stop if nodes < 0)
e = [0.0,0.0,0.0]
for g in range(3):
	#open output file that will contain the wavefunctions  
	f = open("wavefn"+str(g)+".out", "w");
	
	y = [0.0 for i in range(mesh+1)]
	e[g] = solve_sheq(n[g], l, zeta, mesh, dx, r, sqr, r2, vpot, y)
	f.write(str(e[g])+"\n")
	for i in range(mesh+1):
		f.write(str(y[i])+" ")
	f.write("\n")

#####################################################################################################################

print "\n\nzeroth energy of Lithum in eV"
print "1s: ",e[0]*13.6
print "1s:",e[1]*13.6
print "2s:",e[2]*13.6

#####################################################################################################################
#yy = [y,y,y]
w=1
it=0
while True:
	ei=[0.0,0.0,0.0]
	#print "\n\n"
	p = open("wavefn0.out");
	q = open("wavefn1.out");
	s = open("wavefn2.out");

	with p as f:
		inp = f.readlines()

	inp = [x.strip().split(" ") for x in inp]	
	e0 = float(inp[0][0])
	R0=map(float, inp[1])
	p.close()

	with q as f:
		inp = f.readlines()

	inp = [x.strip().split(" ") for x in inp]	
	e1 = float(inp[0][0])
	R1=map(float, inp[1])
	q.close()

	with s as f:
		inp = f.readlines()

	inp = [x.strip().split(" ") for x in inp]	
	e2 = float(inp[0][0])
	R2=map(float, inp[1])
	#print "b",len(R2)
	s.close()
		
	ef = [e0,e1,e2]
	R = [R0,R1,R2]
	
	vhx=[0 for i in range(mesh+1)]
	for i in xrange(3):
		t=0	
		ppp=copy.deepcopy(vpot)
		for j in range(3):
			if j<i:
				
				vhx=v_of_rho(mesh,dx,r,r2,R[j])
				for m in range(mesh+1):
					ppp[m]=ppp[m]+vhx[m]
					
		yy=[0 for x in range(mesh+1)]				
		ei[i] = solve_sheq(n[i], l, zeta, mesh, dx, r, sqr, r2, ppp, yy)
		
		b = open("wavefn"+str(i)+".out", "w");
	
		b.write(str(ei[i])+"\n")
		for j in range(mesh+1):
			b.write(str(yy[j])+" ")
		b.write("\n")
		b.close()
	if(abs(ei[0]-ef[0])<1e-6):
		t+=1
	if(abs(ei[1]-ef[1])<1e-6):
		t+=1
	if(abs(ei[2]-ef[2])<1e-6):
		t+=1	
	
	if(t==3):
		break
	w+=1


#output
####################################################################################
#ionisation evergy
print "\n\nionisation evergy of Lithum in Ry Units"
print "1st IE: ",ei[0]
print "2nd IE:",ei[1]
print "3rd IE:",ei[2]

#ionisation evergy
print "\n\nionisation evergy of Lithum in eV"
print "1st IE: ",ei[0]*13.6
print "2nd IE:",ei[1]*13.6
print "3rd IE:",ei[2]*13.6

####################################################################################	
print "\n\nerror in ionisation energy actual and derived"
print "1s electron:",abs(-9-ei[0])/9*100,"%"
print "1s 2nd electron:",abs(-5.56177-ei[1])/5.56177*100,"%"
print "2s electron:",abs(-0.39645-ei[2])/0.39645*100,"%"
#0.39645 	5.56177  	9

#####################################################################################################################
# solving for ground state lithium wihtout removing the electrons
#yy = [y,y,y]
w=1
it=[y,y,y]
while True:
	ei=[0.0,0.0,0.0]
	#print "\n\n"
	p = open("wavefn0.out");
	q = open("wavefn1.out");
	s = open("wavefn2.out");

	with p as f:
		inp = f.readlines()

	inp = [x.strip().split(" ") for x in inp]	
	e0 = float(inp[0][0])
	R0=map(float, inp[1])
	p.close()

	with q as f:
		inp = f.readlines()

	inp = [x.strip().split(" ") for x in inp]	
	e1 = float(inp[0][0])
	R1=map(float, inp[1])
	q.close()

	with s as f:
		inp = f.readlines()

	inp = [x.strip().split(" ") for x in inp]	
	e2 = float(inp[0][0])
	R2=map(float, inp[1])
	#print "b",len(R2)
	s.close()
		
	ef = [e0,e1,e2]
	R = [R0,R1,R2]
	
	vhx=[0 for i in range(mesh+1)]
	for i in xrange(3):
		t=0	
		ppp=copy.deepcopy(vpot)
		for j in range(3):
			if j!=i:
				
				vhx=v_of_rho(mesh,dx,r,r2,R[j])
				for m in range(mesh+1):
					ppp[m]=ppp[m]+vhx[m]
					
		yy=[0 for x in range(mesh+1)]				
		ei[i] = solve_sheq(n[i], l, zeta, mesh, dx, r, sqr, r2, ppp, yy)
		it[i] = yy
		
		b = open("wavefn"+str(i)+".out", "w");
	
		b.write(str(ei[i])+"\n")
		for j in range(mesh+1):
			b.write(str(yy[j])+" ")
		b.write("\n")
		b.close()
	if(abs(ei[0]-ef[0])<1e-6):
		t+=1
	if(abs(ei[1]-ef[1])<1e-6):
		t+=1
	if(abs(ei[2]-ef[2])<1e-6):
		t+=1	
	
	if(t==3):
		'''plt.plot(r,yy)
		plt.show()'''
		break
	w+=1

#output
####################################################################################
#evergy eigenvalue
print "\n\nevergy eigenvalue of Lithum in Ry Units"
print "1s electron energy: ",ei[0]
print "1s 2nd electron energy",ei[1]
print "2s electron energy",ei[2]

print "\n\nevergy eigenvalue of Lithum in eV"
print "1s electron energy: ",ei[0]*13.6
print "1s electron energy",ei[1]*13.6
print "2s electron energy",ei[2]*13.6


"""Note: if we take range of k as 2--> get energy value after removal of 2p electron
							k as 1---> get energy value after removal of 2s2p electrons
							hence get the ionisation potencials back ...hence verified"""

############################################################################################

X = [[0.0 for i in xrange(mesh+1)] for j in xrange(3)]
for j in xrange(3):
	for i in xrange(mesh+1):
		X[j][i]=(it[j][i]*sqr[i])

#plots
plt.plot(r,"r")
plt.show()
plt.scatter(vpot,r)
plt.show()
#wave func
plt.plot(r,X[0])
plt.show()
plt.plot(r,X[1])
plt.show()
plt.plot(r,X[2])
plt.show()

