# computer-algebraic determination of energy-momentum conservation [symbolic_energy_mom_conservation.py]
#
# basic coordinates and functions appearing
from sympy import *
i, j, k, l, m = symbols('i j k l m',cls=Idx); x = IndexedBase('x'); eps = symbols('eps')
a = Function('a')(x[0]); H = Function('H')(x[0])
barphi = Function('barphi')(x[0]); bare = Function('bare')(x[0]); barp = Function('barp')(x[0]) 
h1 = Function('h1')(x[0],x[1],x[2],x[3]); h2 = Function('h2')(x[0],x[1],x[2],x[3])
h3 = Function('h3')(x[0],x[1],x[2],x[3])
h0 = Function('h0')(x[0],x[1],x[2],x[3]); hD = Function('hD')(x[0],x[1],x[2],x[3])
deltaphi = Function('deltaphi')(x[0],x[1],x[2],x[3])
deltae = Function('deltae')(x[0],x[1],x[2],x[3]); deltap = Function('deltap')(x[0],x[1],x[2],x[3]);
v1 = Function('v1')(x[0],x[1],x[2],x[3]); v2 = Function('v2')(x[0],x[1],x[2],x[3])
v3 = Function('v3')(x[0],x[1],x[2],x[3])
theta11 = Function('theta11')(x[0],x[1],x[2],x[3]); theta22 = Function('theta22')(x[0],x[1],x[2],x[3])
theta33 = Function('theta33')(x[0],x[1],x[2],x[3]); theta12 = Function('theta12')(x[0],x[1],x[2],x[3])
theta13 = Function('theta13')(x[0],x[1],x[2],x[3]); theta23 = Function('theta23')(x[0],x[1],x[2],x[3])
pi11 = Function('pi11')(x[0],x[1],x[2],x[3]); pi22 = Function('pi22')(x[0],x[1],x[2],x[3])
pi33 = Function('pi33')(x[0],x[1],x[2],x[3]); pi12 = Function('pi12')(x[0],x[1],x[2],x[3])
pi13 = Function('pi13')(x[0],x[1],x[2],x[3]); pi23 = Function('pi23')(x[0],x[1],x[2],x[3])

# partial derivative with respect to x[i]
def d(h,i): return diff(h,x[i],evaluate=True)
def dd(h,i,j): return diff(h,x[i],x[j],evaluate=True)            # doesn't symmetrize as should :(
def dds(h,i,j): return (dd(h,i,j)+dd(h,j,i))/2                   # symmetrized
def ddd(h,i,j,k): return diff(h,x[i],x[j],x[k],evaluate=True)    # doesn't symmetrize as should :(
def ddds(h,i,j,k): return (ddd(h,i,j,k)+ddd(h,i,k,j))/2          # symmetrized in last two arguments

# restrict to order(eps) for possibly more efficient execution
def trunc1(h): return (h.expand()).subs([(eps**5,0),(eps**4,0),(eps**3,0),(eps**2,0)])

# input metric with both indices down
a00 = (-1-2*eps*h0)
a01 = (-eps*h1);  a02 = (-eps*h2); a03 = (-eps*h3); a10 = a01; a20 = a02; a30 = a03
a11 = (1+2*eps*(theta11-hD)); a22 = (1+2*eps*(theta22-hD)); a33 = (1+2*eps*(theta33-hD)); 
a12 = (2*eps*theta12); a13 = (2*eps*theta13); a23 = (2*eps*theta23); a21 = a12; a31 = a13; a32 = a23
gdown=a*a*Matrix([[a00,a01,a02,a03],[a10,a11,a12,a13],[a20,a21,a22,a23],[a30,a31,a32,a33]])

# input metric with both indices up
b00 = (-1+2*eps*h0)
b01 = (-eps*h1); b02 = (-eps*h2); b03 = (-eps*h3); b10 = b01; b20 = b02; b30 = b03
b11 = (1+2*eps*(hD-theta11)); b22 = (1+2*eps*(hD-theta22)); b33 = (1+2*eps*(hD-theta33));
b12 = (-2*eps*theta12); b13 = (-2*eps*theta13); b23 = (-2*eps*theta23); b21 = b12; b31 = b13; b32 = b23
gup=1/a/a*Matrix([[b00,b01,b02,b03],[b10,b11,b12,b13],[b20,b21,b22,b23],[b30,b31,b32,b33]])

# christoffel symbols
def dgdown(m,n,i): return d(gdown[m,n],i)
def gamma(r,m,n):
    return Rational(1/2)*Sum( gup[r,i]*(dgdown(i,m,n)+dgdown(i,n,m)-dgdown(m,n,i)),(i,0,3)).doit()

# energy-momentum tensor for perfect fluid plus varphi
t00 = d(barphi,0)**2/2 + eps*d(barphi,0)*d(deltaphi,0) + a*a*(bare + eps*(deltae + 2*h0*bare))
t01 = eps*( d(barphi,0)*(d(deltaphi,1) - d(barphi,0)/2*h1) ) + a*a*eps*(bare*(h1-v1)-barp*v1) 
t02 = eps*( d(barphi,0)*(d(deltaphi,2) - d(barphi,0)/2*h2) ) + a*a*eps*(bare*(h2-v2)-barp*v2) 
t03 = eps*( d(barphi,0)*(d(deltaphi,3) - d(barphi,0)/2*h3) ) + a*a*eps*(bare*(h3-v3)-barp*v3) 
t10 = t01; t20 = t02; t30 = t03
t11 = d(barphi,0)*d(barphi,0)/2 + eps*( ((theta11-hD-h0)*d(barphi,0) + d(deltaphi,0))*d(barphi,0)
      ) +  a*a*( barp + eps*(deltap + 2*(theta11-hD)*barp + pi11) )
t22 = d(barphi,0)*d(barphi,0)/2 + eps*( ((theta22-hD-h0)*d(barphi,0) + d(deltaphi,0))*d(barphi,0)
      ) +  a*a*( barp + eps*(deltap + 2*(theta22-hD)*barp + pi22) )
t33 = d(barphi,0)*d(barphi,0)/2 + eps*( ((theta33-hD-h0)*d(barphi,0) + d(deltaphi,0))*d(barphi,0)
      ) +  a*a*( barp + eps*(deltap + 2*(theta33-hD)*barp + pi33) )
t12 = eps*theta12*d(barphi,0)*d(barphi,0) + a*a*eps*(2*theta12*barp + pi12)
t13 = eps*theta13*d(barphi,0)*d(barphi,0) + a*a*eps*(2*theta13*barp + pi13)
t23 = eps*theta23*d(barphi,0)*d(barphi,0) + a*a*eps*(2*theta23*barp + pi23)
t21 = t12; t31 = t13; t32 = t23
tmunu=Matrix([[t00,t01,t02,t03],[t10,t11,t12,t13],[t20,t21,t22,t23],[t30,t31,t32,t33]])

# energy-momentum conservation equation
def preconserve(n,r):                      # second is a dummy variable   
    return ( Sum( gup[m,k]*d(tmunu[m,n],k),(m,0,3),(k,0,3) ).doit() - 
             Sum( gup[m,k]*(tmunu[m,l]*gamma(l,n,k)+
                            tmunu[n,l]*gamma(l,m,k)),(m,0,3),(k,0,3),(l,0,3)).doit() )
conserve=trunc1(Matrix(2,1,preconserve))   # fill only 1st spatial component to save time

# go over to hubble rate
conserve = conserve.subs([(Derivative(a,x[0]),a*H),(Derivative(a*H,x[0]),a*(d(H,0)+H**2))
                          ,(Derivative(a,x[0],2),a*(d(H,0)+H**2))])

# test energy-momentum conservation
# component nu=0
testconserve0 = -(d(bare,0)+3*H*(bare+barp))-a**(-2)*d(barphi,0)*(dd(barphi,0,0)+2*H*d(barphi,0)
    )+eps*( +a**(-2)*(
        + d(barphi,0)*(dd(deltaphi,1,1)+dd(deltaphi,2,2)+dd(deltaphi,3,3)-dd(deltaphi,0,0))
        - (dd(barphi,0,0) + 4*H*d(barphi,0))*d(deltaphi,0)
        + d(barphi,0)**2*(d(h0,0)+3*d(hD,0)-d(h1,1)-d(h2,2)-d(h3,3)
                          -d(theta11,0)-d(theta22,0)-d(theta33,0))
        + 2*d(barphi,0)*(dd(barphi,0,0) + 2*H*d(barphi,0))*h0
        )
    -d(deltae,0)-3*H*(deltae+deltap)
    +(bare+barp)*(3*d(hD,0)-d(v1,1)-d(v2,2)-d(v3,3)
                    -d(theta11,0)-d(theta22,0)-d(theta33,0))
    -H*(pi11+pi22+pi33)
  )
print(expand(conserve[0]-testconserve0))
print("energy-momentum conservation tested: nu=0")

# component nu=1
testconserve1 = eps*(
    -a**(-2)*(dd(barphi,0,0)+2*H*d(barphi,0))*d(deltaphi,1)
    +d(deltap,1)+(bare+barp)*(d(h0,1)+d(v1,0)-d(h1,0)+4*H*(v1-h1))
    +(d(bare,0)+d(barp,0))*(v1-h1)+d(pi11,1)+d(pi12,2)+d(pi13,3)
    #### remnant of non-symmetrized derivatives
    +a**(-2)*d(barphi,0)*(dd(deltaphi,1,0)-dd(deltaphi,0,1))
    )
print(expand(conserve[1]-testconserve1))
print("energy-momentum conservation tested: nu=1")
