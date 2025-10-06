# computer-algebraic determination of the Ricci scalar [symbolic_ricci_scalar.py]
#
# basic coordinates and functions appearing
from sympy import *
i, j, k, l = symbols('i j k l',cls=Idx); x = IndexedBase('x'); eps = symbols('eps')
a = Function('a')(x[0]); H = Function('H')(x[0])
h1 = Function('h1')(x[0],x[1],x[2],x[3]); h2 = Function('h2')(x[0],x[1],x[2],x[3])
h3 = Function('h3')(x[0],x[1],x[2],x[3])
h0 = Function('h0')(x[0],x[1],x[2],x[3]); hD = Function('hD')(x[0],x[1],x[2],x[3])
theta11 = Function('theta11')(x[0],x[1],x[2],x[3]); theta22 = Function('theta22')(x[0],x[1],x[2],x[3])
theta33 = Function('theta33')(x[0],x[1],x[2],x[3]); theta12 = Function('theta12')(x[0],x[1],x[2],x[3])
theta13 = Function('theta13')(x[0],x[1],x[2],x[3]); theta23 = Function('theta23')(x[0],x[1],x[2],x[3])

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
a01 = (-eps*h1);  a02 = (-eps*h2); a03 = (-eps*h3)
a10 = a01; a20 = a02; a30 = a03
a11 = (1+2*eps*(theta11-hD)); a22 = (1+2*eps*(theta22-hD)); a33 = (1+2*eps*(theta33-hD)); 
a12 = (2*eps*theta12); a13 = (2*eps*theta13); a23 = (2*eps*theta23); 
a21 = a12; a31 = a13; a32 = a23
gdown=a*a*Matrix([[a00,a01,a02,a03],[a10,a11,a12,a13],[a20,a21,a22,a23],[a30,a31,a32,a33]])

# input metric with both indices up
b00 = (-1+2*eps*h0)
b01 = (-eps*h1); b02 = (-eps*h2); b03 = (-eps*h3) 
b10 = b01; b20 = b02; b30 = b03
b11 = (1+2*eps*(hD-theta11)); b22 = (1+2*eps*(hD-theta22)); b33 = (1+2*eps*(hD-theta33));
b12 = (-2*eps*theta12); b13 = (-2*eps*theta13); b23 = (-2*eps*theta23) 
b21 = b12; b31 = b13; b32 = b23
gup=1/a/a*Matrix([[b00,b01,b02,b03],[b10,b11,b12,b13],[b20,b21,b22,b23],[b30,b31,b32,b33]])

# check inversion
test=(gdown*gup);  print(trunc1(test));  print("inversion checked")

# christoffel symbols
def dgdown(m,n,i): return d(gdown[m,n],i)
def gamma(r,m,n):
    return Rational(1/2)*Sum( gup[r,i]*(dgdown(i,m,n)+dgdown(i,n,m)-dgdown(m,n,i)),(i,0,3)).doit()

# ricci tensor with both indices down
def prericcidown(m,n):
    return ( Sum( d(gamma(k,m,n),k) - d(gamma(k,m,k),n),(k,0,3) ).doit() +
        Sum( gamma(l,m,n)*gamma(k,l,k) - gamma(l,m,k)*gamma(k,n,l),(k,0,3),(l,0,3) ).doit() )
riccidown=trunc1(Matrix(4,4,prericcidown))

# ricci tensor with indices up and down
ricciupdown = trunc1(gup*riccidown)

# ricci scalar
ricciscalar = ricciupdown.trace()

# go over to hubble rate
ricciscalar = ricciscalar.subs([(Derivative(a,x[0]),a*H),(Derivative(a*H,x[0]),a*(d(H,0)+H**2))
                                ,(Derivative(a,x[0],2),a*(d(H,0)+H**2))])

# adjust expression until the difference subtracts to zero 
testricciscalar = a**(-2)*(6*(d(H,0)+H**2) - 2*eps*(
    +3*H*d(h0,0)+6*(d(H,0)+H**2)*h0
    +dd(h0,1,1)+dd(h0,2,2)+dd(h0,3,3)-2*(dd(hD,1,1)+dd(hD,2,2)+dd(hD,3,3))
    +dd(theta11,1,1)+dd(theta22,1,1)+dd(theta33,1,1)
    +dd(theta11,2,2)+dd(theta22,2,2)+dd(theta33,2,2)
    +dd(theta11,3,3)+dd(theta22,3,3)+dd(theta33,3,3)
    +3*dd(hD,0,0)+9*H*d(hD,0) -(dds(h1,0,1)+dds(h2,0,2)+dds(h3,0,3))-3*H*(d(h1,1)+d(h2,2)+d(h3,3))
    -(dd(theta11,0,0)+dd(theta22,0,0)+dd(theta33,0,0))-3*H*(d(theta11,0)+d(theta22,0)+d(theta33,0))
    -(dd(theta11,1,1)+dd(theta22,2,2)+dd(theta33,3,3))
    -2*(dds(theta12,1,2)+dds(theta13,1,3)+dds(theta23,2,3))   )) 
print(expand(a**2*(ricciscalar - testricciscalar)))
print("ricci scalar checked at 1st order")

