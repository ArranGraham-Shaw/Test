from numpy import *
from pylab import *
from scipy import *
from scipy.constants import *
from scipy.special import *
from scipy.integrate import *
from FWHMcalc import FWHM
from scatsimfourierfunctions import fourierxtok, fourierktox, timeoperator, timeevolve
matplotlib.rc('font', size=31)

# Add some github annotation

# Add some more github annotation

# Even more github annotation

def findy(array,value):     #Define a function to match the randomly generate number to an element of the "dis" array in order to call the corresponding angle
    idx = (np.abs(array-value)).argmin()
    return idx
	
def Pscat(L,x,dx,j,phi):     #Define a function that calculates the probability of a scattering event with a given wavefunction
	a=cos(2*pi*x*sin(L)/j)*phi*conjugate(cos(2*pi*x*sin(L)/j)*phi)*dx
	b=sum(a)/(2*pi)
	return b

def A(position):     #A(x) has been found and is defined in terms of the zeroth order Bessel function of the first kind
	return sqrt(0.5-0.5*jn(0,(4*pi*position)/j))

def B(theta):      #Function to update phi in event of a scattering
	return cos(2*pi*x*sin(theta)/j)	
	
	

dim=1000   #n.b. dimensions of i,l,y,dis and the for loop must be the same
dim2=200
d=1e-6
j=1*d
x=linspace(-d,d,dim)
dx=abs(x[0]-x[1])
k=linspace(-180000000*hbar,180000000*hbar,dim)
dk=abs(k[0]-k[1])
D=x[dim-1]-x[0]
phi=ones(dim)/D
phiprob=phi*conjugate(phi)     #Defined again here so the code will run without the while loop below

l=linspace(-pi,pi,dim2)   #l is the range of possible scattering angles
y=zeros(dim2)


#This segment iterates scattering events and updates phi accordingly	
count=0
while count<200:
	count=count+1
	
	for i in linspace(0,dim2-1,dim2):    #This for loop calculates the scattering probability over the range of angles for the given wavefunction pass
		y[i] = Pscat(l[i],x,dx,j,phi)

	dis=cumsum(y*abs(l[1]-l[2]))   #dis is the cumulative sum of scattering probabilities	
	w=uniform()   #w is a randomly generated float between 0 and 1
	if w>dis[dim2-1]:#w represents a random scattering event with the probability of scattering at a certain angle or a no-scattering event occuring being mapped out by the "dis" array by this if else loop
		scatevent=0
	else:
		scatevent=l[findy(dis,w)]
	
	if scatevent == 0 or abs(scatevent) == pi:
		phi = phi*A(x)/sqrt(sum(phi*A(x)*conjugate(phi*A(x))*dx))
	else:
		phi = phi*B(scatevent)/sqrt(sum(phi*B(scatevent)*conjugate(phi*B(scatevent))*dx))
	print "scatevent", scatevent
	if count % 10 == 0:
		print count

print "Count = ",count
phiprob=phi*conjugate(phi)

phik=fourierxtok(phi,x,dx,k,dim)
phikprob=phik*conjugate(phik)


philoc=ones(dim)
philoc[0:dim/2]=0
philoc[dim/2:dim-1]=sqrt(2)*phi[dim/2:dim-1]
philocprob=philoc*conjugate(philoc)
xprob=philocprob#/sum(philocprob)
x2prob=phiprob/sum(phiprob)


expectx=sum(x*xprob*dx)
#print 'Expected value = ', expectx
expectx2=sum(x**2*xprob*dx)
#print 'Expected value of x^2 = ', expectx
varx=expectx2-expectx**2
#print 'Variance of x = ', varx
print 'Standard deviation of x = ', sqrt(varx)

# philoc=ones(dim)
# philoc[0:dim/2]=0
# philoc[dim/2:dim-1]=sqrt(2)*phi[dim/2:dim-1]
# philocprob=philoc*conjugate(philoc)
kprob=phikprob#/sum(phikprob)
k2prob=phikprob/sum(phikprob)

expectk=sum(k*kprob*dk)
#print 'Expected value = ', expectx
expectk2=sum(k**2*kprob*dk)
#print 'Expected value of x^2 = ', expectx
vark=expectk2-expectk**2
#print 'Variance of k = ', vark
print 'Standard deviation of k = ', sqrt(vark)


#This segment calculates the area under the graph to check that phi is normalised
phiprob=phi*conjugate(phi)
area=trapz(phiprob, dx=dx)
print 'Normalisation check of position = ', area
	
phikprob=phik*conjugate(phik)
area=trapz(phikprob, dx=dk)
print 'Normalisation check of momentum = ', area


subplot(2,1,1)
plot(x,x2prob)
title('Relative Position')
xlabel('x[m]')
ylabel('P(x)')
locator_params(nbins=4)
ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
subplot(2,1,2)
plot(k,k2prob)
title('Relative Momentum')
xlabel('p[kg m/s]')
ylabel('P(p)')
locator_params(nbins=4)
ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ticklabel_format(style='sci', axis='y', scilimits=(0,0))
show()
# plot(l,dis)
# show()
