
(*source term for scalar-induced gravitational waves*)

 Clear[h0];
 Clear[hD];
 Clear[h];

(*metric with both indices down*)

 gdown[0,0]:=a[x[0]]^2(-1-2eps h0[x[0],x[1],x[2],x[3]]);
 gdown[0,1]:=a[x[0]]^2(eps Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]);
 gdown[0,2]:=a[x[0]]^2(eps Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]);
 gdown[0,3]:=a[x[0]]^2(eps Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]]);
 gdown[1,0]:=gdown[0,1];
 gdown[2,0]:=gdown[0,2];
 gdown[3,0]:=gdown[0,3];
 gdown[1,1]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]);
 gdown[2,2]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]);
 gdown[3,3]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]);
 gdown[1,2]:=0;
 gdown[1,3]:=0;
 gdown[2,3]:=0;
 gdown[2,1]:=gdown[1,2];
 gdown[3,1]:=gdown[1,3];
 gdown[3,2]:=gdown[2,3];

(*inverse metric with both indices up*)

 gup[0,0]:=a[x[0]]^(-2)(-1+2eps h0[x[0],x[1],x[2],x[3]] + eps^2(
  -4 h0[x[0],x[1],x[2],x[3]]^2 
  +  Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]^2
  +  Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]^2
  +  Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]]^2
 ));
 gup[0,1]:=a[x[0]]^(-2) Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]] (eps + 
 eps^2(
  2 hD[x[0],x[1],x[2],x[3]] 
 -2 h0[x[0],x[1],x[2],x[3]] 
 ));
 gup[0,2]:=a[x[0]]^(-2) Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]] (eps + 
 eps^2(
  2 hD[x[0],x[1],x[2],x[3]] 
 -2 h0[x[0],x[1],x[2],x[3]] 
 ));
 gup[0,3]:=a[x[0]]^(-2) Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]] (eps + 
 eps^2(
  2 hD[x[0],x[1],x[2],x[3]] 
 -2 h0[x[0],x[1],x[2],x[3]] 
 ));
 gup[1,0]:=gup[0,1];
 gup[2,0]:=gup[0,2];
 gup[3,0]:=gup[0,3];
 gup[1,1]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]] +eps^2(
   4 hD[x[0],x[1],x[2],x[3]]^2 
 - Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]] 
   Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]
 ));
 gup[2,2]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]] +eps^2(
   4 hD[x[0],x[1],x[2],x[3]]^2 
 - Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]] 
   Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]
 ));
 gup[3,3]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]] +eps^2(
   4 hD[x[0],x[1],x[2],x[3]]^2 
 - Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]] 
   Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]]
 ));
 gup[1,2]:=a[x[0]]^(-2)(0+eps^2(
 - Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]] 
   Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]
 ));
 gup[1,3]:=a[x[0]]^(-2)(0+eps^2(
 - Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]] 
   Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]]
 ));
 gup[2,3]:=a[x[0]]^(-2)(0+eps^2(
 - Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]] 
   Derivative[0,0,0,1][h][x[0],x[1],x[2],x[3]]
 ));
 gup[2,1]:=gup[1,2];
 gup[3,1]:=gup[1,3];
 gup[3,2]:=gup[2,3];

(*check inversion*)

  product[mu_,nu_]:=Sum[gdown[mu,sigma] gup[sigma,nu],{sigma,0,3}]
  Print[Series[product[0,0],{eps,0,2}]]
  Print[Series[product[1,1],{eps,0,2}]]
  Print[Series[product[2,2],{eps,0,2}]]
  Print[Series[product[3,3],{eps,0,2}]]
  Print[Series[product[0,1],{eps,0,2}]]
  Print[Series[product[0,2],{eps,0,2}]]
  Print[Series[product[0,3],{eps,0,2}]]
  Print[Series[product[1,2],{eps,0,2}]]
  Print[Series[product[1,3],{eps,0,2}]]
  Print[Series[product[2,3],{eps,0,2}]]
  Print["inversion checked"]

(*cristoffel symbols*)

 gamma[rho_,mu_,nu_]:=1/2 Sum[gup[rho,sigma](D[gdown[sigma,mu],x[nu]]+
  D[gdown[sigma,nu],x[mu]]-D[gdown[mu,nu],x[sigma]]),{sigma,0,3}]

(*ricci tensor with both indices down*)

 riccidown[mu_,nu_]:=Sum[D[gamma[alpha,mu,nu],x[alpha]],{alpha,0,3}]-
  Sum[D[gamma[alpha,mu,alpha],x[nu]],{alpha,0,3}]+
  Sum[gamma[beta,mu,nu]gamma[alpha,beta,alpha],{alpha,0,3},{beta,0,3}]-
  Sum[gamma[beta,mu,alpha]gamma[alpha,beta,nu],{alpha,0,3},{beta,0,3}]

(*raise first index in ricci tensor*)

 ricciupdown[mu_,nu_]:=Sum[gup[mu,alpha]riccidown[alpha,nu],{alpha,0,3}]

(*ricci scalar*)

 ricciscalar:=Sum[ricciupdown[alpha,alpha],{alpha,0,3}]

(*einstein tensor with both indices down*)

 einstein[mu_,nu_]:=riccidown[mu,nu] - ricciscalar/2 gdown[mu,nu]  

(*adjust expression until the difference subtracts to zero*)

 testeinstein[1,2]:=(0 + eps(
 - Derivative[0,1,1,0][h0][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,1,0][hD][x[0],x[1],x[2],x[3]]
 - Derivative[1,1,1,0][h][x[0],x[1],x[2],x[3]]
 - 2 H[x[0]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 )
		     + eps^2(
 + 2 h0[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h0][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,0,0][h0][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h0][x[0],x[1],x[2],x[3]] 
 + 2 hD[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][hD][x[0],x[1],x[2],x[3]]
 + 3 Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]] 
 - Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h0][x[0],x[1],x[2],x[3]] 
 - Derivative[0,1,0,0][h0][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]] 
 + Derivative[1,0,0,0][h0][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 + 4 H[x[0]] h0[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 + 2 h0[x[0],x[1],x[2],x[3]] Derivative[1,1,1,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[1,0,1,0][h][x[0],x[1],x[2],x[3]] 
 - Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]] Derivative[1,1,0,0][h][x[0],x[1],x[2],x[3]]
 - 2 H[x[0]] Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]
 - 2 H[x[0]] Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[1,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[1,0,1,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]] Derivative[0,0,0,2][h][x[0],x[1],x[2],x[3]]
 - Derivative[0,1,0,1][h][x[0],x[1],x[2],x[3]] Derivative[0,0,1,1][h][x[0],x[1],x[2],x[3]] 
 )
 ); 

 subst={a''[x[0]]->a[x[0]](H'[x[0]]+H[x[0]]^2),a'[x[0]]->a[x[0]] H[x[0]]};
 Print[Expand[Normal[Series[Expand[(einstein[1,2]//.subst)]-testeinstein[1,2],{eps,0,2}]]]]
 Print["einstein 12 checked to 2nd order"]

(*construct effective energy-momentum tensor by moving metric perturbations to right-hand side*)

 rhs[1,2]:=-testeinstein[1,2]+eps^2(
 -2 Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h0][x[0],x[1],x[2],x[3]] 
 -2 Derivative[0,1,0,0][h0][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]] 
 +2 (H'[x[0]]/H[x[0]]^2 - 1) (
    Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]])
 -2/H[x[0]] (
    Derivative[1,1,0,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][hD][x[0],x[1],x[2],x[3]]
  + Derivative[1,0,1,0][hD][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][hD][x[0],x[1],x[2],x[3]]
 )
 );

(*substitute bardeen potentials*)

 h0[x0_,x1_,x2_,x3_]:=(
   phi[x0,x1,x2,x3]
 - Derivative[1,0,0,0][h][x0,x1,x2,x3]
 - H[x0] h[x0,x1,x2,x3]);
 
 hD[x0_,x1_,x2_,x3_]:=(
   psi[x0,x1,x2,x3]
 + H[x0] h[x0,x1,x2,x3]);

(*adjust expressions until the difference subtracts to zero*)

 testrhs[1,2]:=(0 + eps (
   Derivative[0,1,1,0][phi][x[0],x[1],x[2],x[3]]
 - Derivative[0,1,1,0][psi][x[0],x[1],x[2],x[3]]
     )
 + eps^2(
 - Derivative[0,1,0,0][phi][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][phi][x[0],x[1],x[2],x[3]]  
 - 2 phi[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][phi][x[0],x[1],x[2],x[3]] 
 - Derivative[0,1,0,0][phi][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][psi][x[0],x[1],x[2],x[3]]  
 - Derivative[0,0,1,0][phi][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][psi][x[0],x[1],x[2],x[3]]  
 - 2 psi[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][psi][x[0],x[1],x[2],x[3]] 
 + (2 H'[x[0]]/H[x[0]]^2 - 5 )( 
   Derivative[0,1,0,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][psi][x[0],x[1],x[2],x[3]] ) 
 - 2/H[x[0]]( 
  +Derivative[1,1,0,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][psi][x[0],x[1],x[2],x[3]] 
  +Derivative[1,0,1,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][psi][x[0],x[1],x[2],x[3]] 
   ) 
 (* *)
 + 2 H[x[0]](
  +Derivative[0,1,1,0][phi][x[0],x[1],x[2],x[3]] h[x[0],x[1],x[2],x[3]]
  -phi[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 )
 - Derivative[1,0,0,0][phi][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 + 2 Derivative[0,1,1,0][phi][x[0],x[1],x[2],x[3]] Derivative[1,0,0,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,0,0][phi][x[0],x[1],x[2],x[3]] Derivative[1,0,1,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,1,0][phi][x[0],x[1],x[2],x[3]] Derivative[1,1,0,0][h][x[0],x[1],x[2],x[3]]
 (* *)
 -2 H[x[0]](
  +Derivative[0,1,0,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]
  +Derivative[0,0,1,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]
  +psi[x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
  +Derivative[0,1,1,0][psi][x[0],x[1],x[2],x[3]] h[x[0],x[1],x[2],x[3]]
 )
 - Derivative[1,0,0,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[1,0,1,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,1,0,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[1,1,0,0][psi][x[0],x[1],x[2],x[3]] Derivative[0,0,1,0][h][x[0],x[1],x[2],x[3]]
 (* *)
 + Derivative[0,2,0,0][h][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]] Derivative[0,0,2,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,0,1][h][x[0],x[1],x[2],x[3]] Derivative[0,0,1,1][h][x[0],x[1],x[2],x[3]]
 - Derivative[0,2,0,0][h][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[0,0,2,0][h][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[0,0,0,2][h][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 + Derivative[2,0,0,0][h][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
 - Derivative[1,0,1,0][h][x[0],x[1],x[2],x[3]] Derivative[1,1,0,0][h][x[0],x[1],x[2],x[3]]
 + 2 H[x[0]](
  +Derivative[1,0,0,0][h][x[0],x[1],x[2],x[3]] Derivative[0,1,1,0][h][x[0],x[1],x[2],x[3]]
  )
 )
 );

 Print[Expand[Normal[Series[Expand[rhs[1,2]]-testrhs[1,2],{eps,0,2}]]]]
 Print["sigw 12 checked to 2nd order"]
