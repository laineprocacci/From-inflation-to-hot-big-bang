# computer-algebraic determination of the source for scalar-induced GWs [symbolic_sigw.py]
#
# basic coordinates and functions appearing
from sympy import *
i, j, k, l = symbols('i j k l',cls=Idx); x = IndexedBase('x'); eps = symbols('eps')
a = Function('a')(x[0]); H = Function('H')(x[0])
h = Function('h')(x[0],x[1],x[2],x[3])
h0 = Function('h0')(x[0],x[1],x[2],x[3]); hD = Function('hD')(x[0],x[1],x[2],x[3])
phi = Function('phi')(x[0],x[1],x[2],x[3]); psi = Function('psi')(x[0],x[1],x[2],x[3])

# partial derivative with respect to x[i]
def d(h,i): return diff(h,x[i],evaluate=True)
def dd(h,i,j): return diff(h,x[i],x[j],evaluate=True)            # doesn't symmetrize as should :(
def dds(h,i,j): return (dd(h,i,j)+dd(h,j,i))/2                   # symmetrized
def ddd(h,i,j,k): return diff(h,x[i],x[j],x[k],evaluate=True)    # doesn't symmetrize as should :(
def ddds(h,i,j,k): return (ddd(h,i,j,k)+ddd(h,i,k,j))/2          # symmetrized in last two args

# restrict to eps^order for possibly more efficient execution
def trunc2(h): return (h.expand()).subs([(eps**5,0),(eps**4,0),(eps**3,0)])
def trunc1(h): return (h.expand()).subs([(eps**2,0)])
order = 2

# input metric with both indices down
a00 = a**2*(-1-2*eps*h0)
a01 = a**2*(eps*d(h,1));  a02 = a**2*(eps*d(h,2)); a03 = a**2*(eps*d(h,3))
a10 = a01; a20 = a02; a30 = a03
a11 = a**2*(1-2*eps*hD); a22=a11; a33=a11
a12 = 0; a13 = 0; a23 = 0; a21 = a12; a31 = a13; a32 = a23
gdown=Matrix([[a00,a01,a02,a03],[a10,a11,a12,a13],[a20,a21,a22,a23],[a30,a31,a32,a33]])

# input metric with both indices up
b00 = a**(-2)*(-1+2*eps*h0 + eps**2*( -4*h0**2 + d(h,1)**2 + d(h,2)**2 + d(h,3)**2 ))
b01 = a**(-2)*d(h,1)*(eps +  eps**2*(  2*hD - 2*h0  ))
b02 = a**(-2)*d(h,2)*(eps +  eps**2*(  2*hD - 2*h0  ))
b03 = a**(-2)*d(h,3)*(eps +  eps**2*(  2*hD - 2*h0  ))
b10 = b01; b20 = b02; b30 = b03
b11 = a**(-2)*(1 + 2*eps*hD + eps**2*( 4*hD**2 - d(h,1)**2 ))
b22 = a**(-2)*(1 + 2*eps*hD + eps**2*( 4*hD**2 - d(h,2)**2 ))
b33 = a**(-2)*(1 + 2*eps*hD + eps**2*( 4*hD**2 - d(h,3)**2 ))
b12 = a**(-2)*(0 + eps**2*( - d(h,1)*d(h,2) ))
b13 = a**(-2)*(0 + eps**2*( - d(h,1)*d(h,3) ))
b23 = a**(-2)*(0 + eps**2*( - d(h,2)*d(h,3) ))
b21 = b12; b31 = b13; b32 = b23
gup=Matrix([[b00,b01,b02,b03],[b10,b11,b12,b13],[b20,b21,b22,b23],[b30,b31,b32,b33]])
if order==1: gup = trunc1(gup) 

# check inversion
test=(gdown*gup) 
if order==1: print(trunc1(test)) 
else: print(trunc2(test)) 
print("inversion checked")

# christoffel symbols
def dgdown(m,n,i): return d(gdown[m,n],i)
def gamma(r,m,n):
    return Rational(1/2)*Sum( gup[r,i]*(dgdown(i,m,n)+dgdown(i,n,m)-dgdown(m,n,i)),(i,0,3)).doit()

# ricci tensor with both indices down
def prericcidown(m,n):
    return ( Sum( d(gamma(k,m,n),k) - d(gamma(k,m,k),n),(k,0,3) ).doit() +
        Sum( gamma(l,m,n)*gamma(k,l,k) - gamma(l,m,k)*gamma(k,l,n),(k,0,3),(l,0,3) ).doit() )
riccidown=trunc2(Matrix(4,4,prericcidown))
if order==1: riccidown=trunc1(riccidown)

# ricci tensor with indices up and down
ricciupdown = trunc2(gup*riccidown)
if order==1: ricciupdown=trunc1(ricciupdown)

# ricci scalar
ricciscalar = ricciupdown.trace()

# einstein tensor with both indices down
einstein = trunc2(riccidown - ricciscalar*gdown/2)
if order==1: einstein=trunc1(einstein)
    
# replace derivatives of the scale factor by the hubble rate
einstein = einstein.subs([(Derivative(a,x[0]),a*H),(Derivative(a,x[0],2),a*(d(H,0)+H**2))])

# adjust expression until the difference subtracts to zero 
testeinstein12 = (0 + eps*( - dds(h0,1,2) + dds(hD,1,2) - ddds(h,0,1,2) - 2*H*dds(h,1,2) )
                  + eps**2*(
                      + d(h0,1)*d(h0,2) + 3*d(hD,1)*d(hD,2)
                      + 2*( h0*dds(h0,1,2) + hD*dds(hD,1,2) )
                      - d(h0,1)*d(hD,2) - d(hD,1)*d(h0,2)
                      + ( d(h0,0) + d(hD,0) )*dds(h,1,2)
                      + 4*H*h0*dds(h,1,2) + 2*h0*ddds(h,0,1,2)
                      - d(hD,1)*dd(h,0,2) - d(hD,2)*dd(h,0,1)
                      - dd(hD,0,1)*d(h,2) - dd(hD,0,2)*d(h,1)
                      - 2*H*( d(hD,1)*d(h,2) + d(hD,2)*d(h,1) )
                      + ( dd(h,1,1) + dd(h,2,2) + dd(h,3,3))*dds(h,1,2)
                      - ( dd(h,1,1)*dds(h,1,2) + dds(h,2,1)*dd(h,2,2) + dd(h,3,1)*dd(h,3,2) )
    #### remnant of non-symmetrized derivatives
                      + (dd(h,1,1)-dd(h,2,2))*(dd(h,2,1)-dd(h,1,2))/2   ) ) 
if order==1: testeinstein12=trunc1(testeinstein12)
print(expand((einstein[1,2]+einstein[2,1])/2 - testeinstein12))
print("einstein 12 checked to 2nd order")

# construct effective energy-momentum tensor by moving metric perturbations to right-hand side
rhs12=-testeinstein12 + eps**2*(
    -2*( d(h0,1)*d(hD,2) + d(hD,1)*d(h0,2) ) 
    -2/H*( dd(hD,0,1)*d(hD,2) + dd(hD,0,2)*d(hD,1) )
    +2*(d(H,0)/H**2 - 1)*d(hD,1)*d(hD,2) )
if order==1: rhs12=trunc1(rhs12)

# substitute bardeen potentials
rhs12=rhs12.subs([(h0,phi-d(h,0)-H*h),(hD,psi+H*h)]).doit()
if order==1: rhs12=trunc1(rhs12)

# adjust expressions until the difference subtracts to zero
testrhs12 = (0 + eps*( dds(phi,1,2)-dds(psi,1,2) )
             + eps**2*(
                 -d(phi,1)*d(phi,2) - 2*phi*dds(phi,1,2)
                 -( d(phi,1)*d(psi,2) + d(phi,2)*d(psi,1) ) - 2*psi*dds(psi,1,2)
                 +( 2*d(H,0)/H**2 - 5 )*d(psi,1)*d(psi,2)
                 -2/H*( dd(psi,0,1)*d(psi,2) + dd(psi,0,2)*d(psi,1) )
                 +2*H*( h*dds(phi,1,2) - phi*dds(h,1,2) )
                 -d(phi,0)*dds(h,1,2) + 2*dds(phi,1,2)*d(h,0) + d(phi,1)*dd(h,0,2) + d(phi,2)*dd(h,0,1)
                 -2*H*( d(psi,1)*d(h,2) + d(psi,2)*d(h,1) + psi*dds(h,1,2) + h*dds(psi,1,2) )
                 -d(psi,0)*dds(h,1,2) - dd(psi,0,1)*d(h,2) - dd(psi,0,2)*d(h,1)
                 +dd(h,1,1)*dds(h,1,2) + dds(h,2,1)*dd(h,2,2) + dd(h,3,1)*dd(h,3,2)
                 -( dd(h,1,1) + dd(h,2,2) + dd(h,3,3))*dds(h,1,2)
                 +dd(h,0,0)*dds(h,1,2) - dd(h,0,1)*dd(h,0,2) + 2*H*d(h,0)*dds(h,1,2)
    #### remnant of non-symmetrized derivatives
                 - (dd(h,1,1)-dd(h,2,2))*(dd(h,2,1)-dd(h,1,2))/2  ) )
if order==1: testrhs12=trunc1(testrhs12)
print(expand(rhs12-testrhs12))
print("sigw 12 checked to 2nd order")
