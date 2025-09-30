
(*ricci scalar to first order in the presence of a general perturbation*)

(*metric*)

 gdown[0,0]:=a[x[0]]^2(-1-2eps h0[x[0],x[1],x[2],x[3]]);
 gdown[0,1]:=a[x[0]]^2(-eps h1[x[0],x[1],x[2],x[3]]);
 gdown[0,2]:=a[x[0]]^2(-eps h2[x[0],x[1],x[2],x[3]]);
 gdown[0,3]:=a[x[0]]^2(-eps h3[x[0],x[1],x[2],x[3]]);
 gdown[1,0]:=gdown[0,1];
 gdown[2,0]:=gdown[0,2];
 gdown[3,0]:=gdown[0,3];
 gdown[1,1]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]+2 eps theta11[x[0],x[1],x[2],x[3]]);
 gdown[2,2]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]+2 eps theta22[x[0],x[1],x[2],x[3]]);
 gdown[3,3]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]+2 eps theta33[x[0],x[1],x[2],x[3]]);
 gdown[1,2]:=a[x[0]]^2(2eps theta12[x[0],x[1],x[2],x[3]]);
 gdown[1,3]:=a[x[0]]^2(2eps theta13[x[0],x[1],x[2],x[3]]);
 gdown[2,3]:=a[x[0]]^2(2eps theta23[x[0],x[1],x[2],x[3]]);
 gdown[2,1]:=gdown[1,2];
 gdown[3,1]:=gdown[1,3];
 gdown[3,2]:=gdown[2,3];

(*inverse metric*)

 gup[0,0]:=a[x[0]]^(-2)(-1+2eps h0[x[0],x[1],x[2],x[3]]);
 gup[0,1]:=a[x[0]]^(-2)(-eps h1[x[0],x[1],x[2],x[3]]);
 gup[0,2]:=a[x[0]]^(-2)(-eps h2[x[0],x[1],x[2],x[3]]);
 gup[0,3]:=a[x[0]]^(-2)(-eps h3[x[0],x[1],x[2],x[3]]);
 gup[1,0]:=gup[0,1];
 gup[2,0]:=gup[0,2];
 gup[3,0]:=gup[0,3];
 gup[1,1]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]]-2 eps theta11[x[0],x[1],x[2],x[3]]);
 gup[2,2]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]]-2 eps theta22[x[0],x[1],x[2],x[3]]);
 gup[3,3]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]]-2 eps theta33[x[0],x[1],x[2],x[3]]);
 gup[1,2]:=a[x[0]]^(-2)(-2eps theta12[x[0],x[1],x[2],x[3]]);
 gup[1,3]:=a[x[0]]^(-2)(-2eps theta13[x[0],x[1],x[2],x[3]]);
 gup[2,3]:=a[x[0]]^(-2)(-2eps theta23[x[0],x[1],x[2],x[3]]);
 gup[2,1]:=gup[1,2];
 gup[3,1]:=gup[1,3];
 gup[3,2]:=gup[2,3];

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

(*adjust expression until the difference subtracts to zero*)

  calH=a'[x[0]]/a[x[0]];
  calHp=a''[x[0]]/a[x[0]] - (a'[x[0]]/a[x[0]])^2;
  testricciscalar=6 a''[x[0]]/a[x[0]]^3-2 eps/a[x[0]]^2(  
      + 3 calH Derivative[1,0,0,0][h0][x[0],x[1],x[2],x[3]]
      + 6 (calHp+calH^2) Derivative[0,0,0,0][h0][x[0],x[1],x[2],x[3]]
      + Derivative[0,2,0,0][h0][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][h0][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][h0][x[0],x[1],x[2],x[3]]
      - 2 Derivative[0,2,0,0][hD][x[0],x[1],x[2],x[3]]
      - 2 Derivative[0,0,2,0][hD][x[0],x[1],x[2],x[3]]
      - 2 Derivative[0,0,0,2][hD][x[0],x[1],x[2],x[3]]
      + Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]
      + Derivative[0,2,0,0][theta22][x[0],x[1],x[2],x[3]]
      + Derivative[0,2,0,0][theta33][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][theta11][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][theta22][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][theta33][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][theta11][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][theta22][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][theta33][x[0],x[1],x[2],x[3]]
      + 3 Derivative[2,0,0,0][hD][x[0],x[1],x[2],x[3]]  
      + 9 calH Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]  
      - Derivative[1,1,0,0][h1][x[0],x[1],x[2],x[3]]
      - Derivative[1,0,1,0][h2][x[0],x[1],x[2],x[3]]
      - Derivative[1,0,0,1][h3][x[0],x[1],x[2],x[3]]
      - 3 calH Derivative[0,1,0,0][h1][x[0],x[1],x[2],x[3]]
      - 3 calH Derivative[0,0,1,0][h2][x[0],x[1],x[2],x[3]]
      - 3 calH Derivative[0,0,0,1][h3][x[0],x[1],x[2],x[3]]
      - Derivative[2,0,0,0][theta11][x[0],x[1],x[2],x[3]]
      - Derivative[2,0,0,0][theta22][x[0],x[1],x[2],x[3]]
      - Derivative[2,0,0,0][theta33][x[0],x[1],x[2],x[3]]
      - 3 calH Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
      - 3 calH Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
      - 3 calH Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
      - Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]
      - Derivative[0,0,2,0][theta22][x[0],x[1],x[2],x[3]]
      - Derivative[0,0,0,2][theta33][x[0],x[1],x[2],x[3]]
      - 2 Derivative[0,1,1,0][theta12][x[0],x[1],x[2],x[3]]
      - 2 Derivative[0,1,0,1][theta13][x[0],x[1],x[2],x[3]]
      - 2 Derivative[0,0,1,1][theta23][x[0],x[1],x[2],x[3]]
    );

 Print[Series[ricciscalar-testricciscalar,{eps,0,1}]]
 Print["ricci scalar checked at 1st order"]
