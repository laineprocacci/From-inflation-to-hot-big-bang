
(*einstein tensor to first order in the presence of a general perturbation*)

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

(*print out zeroth order*)
(*
 Do[Print[rho,mu,nu,": ",Series[gamma[rho,mu,nu],{eps,0,0}]," "],{rho,0,3},{mu,0,3},{nu,0,3}]
 Abort[]
 *)

(*ricci tensor*)

 riccidown[mu_,nu_]:=Sum[D[gamma[alpha,mu,nu],x[alpha]],{alpha,0,3}]-
  Sum[D[gamma[alpha,mu,alpha],x[nu]],{alpha,0,3}]+
  Sum[gamma[beta,mu,nu]gamma[alpha,beta,alpha],{alpha,0,3},{beta,0,3}]-
  Sum[gamma[beta,mu,alpha]gamma[alpha,nu,beta],{alpha,0,3},{beta,0,3}]

(*print out zeroth order*)
(*
 Do[Print[mu,nu,": ",Series[riccidown[mu,nu],{eps,0,0}]," "],{mu,0,3},{nu,0,3}]
 Abort[]
 *)

(*raise first index in ricci tensor*)

 ricciupdown[mu_,nu_]:=Sum[gup[mu,alpha]riccidown[alpha,nu],{alpha,0,3}]

(*ricci scalar*)

 ricciscalar:=Sum[ricciupdown[alpha,alpha],{alpha,0,3}]

(*print out zeroth order*)
(*
 Do[Print["R: ",Series[ricciscalar,{eps,0,0}]," "],{mu,0,0},{nu,0,0}]
 Abort[]
 *)

(*einstein tensor with both indices down*)

 einstein[mu_,nu_]:=riccidown[mu,nu] - ricciscalar/2 gdown[mu,nu]  

(*adjust expression until the difference subtracts to zero*)

  H=a'[x[0]]/a[x[0]];
  Hp=a''[x[0]]/a[x[0]] - (a'[x[0]]/a[x[0]])^2;

  testeinstein[0,0]:=(3 H^2 + eps(
   2(
      + Derivative[0,2,0,0][hD][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][hD][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][hD][x[0],x[1],x[2],x[3]]
      - 3 H Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
      + H (  Derivative[0,1,0,0][h1][x[0],x[1],x[2],x[3]]
           + Derivative[0,0,1,0][h2][x[0],x[1],x[2],x[3]]
           + Derivative[0,0,0,1][h3][x[0],x[1],x[2],x[3]] 
	     ))
      + Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][theta22][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][theta33][x[0],x[1],x[2],x[3]]
      + 2 Derivative[0,1,1,0][theta12][x[0],x[1],x[2],x[3]]
      + 2 Derivative[0,1,0,1][theta13][x[0],x[1],x[2],x[3]]
      + 2 Derivative[0,0,1,1][theta23][x[0],x[1],x[2],x[3]]
      - (
        + Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,2,0,0][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[0,2,0,0][theta33][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta33][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta33][x[0],x[1],x[2],x[3]]
        )
      + 2 H ( Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
            + Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
	    + Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]] 
            )
   ))

  Print["component 00:"]
  Print[Series[einstein[0,0]-testeinstein[0,0],{eps,0,1}]]
  (* Abort[] *)

  testeinstein[0,1]:=(0 + eps(
   2(
      + H Derivative[0,1,0,0][h0][x[0],x[1],x[2],x[3]]
      + Derivative[1,1,0,0][hD][x[0],x[1],x[2],x[3]]
      )
   +(2 Hp + H^2) Derivative[0,0,0,0][h1][x[0],x[1],x[2],x[3]]
   +1/2(
      + Derivative[0,2,0,0][h1][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][h1][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][h1][x[0],x[1],x[2],x[3]]
      - Derivative[0,2,0,0][h1][x[0],x[1],x[2],x[3]]
      - Derivative[0,1,1,0][h2][x[0],x[1],x[2],x[3]]
      - Derivative[0,1,0,1][h3][x[0],x[1],x[2],x[3]]
      )
      + Derivative[1,1,0,0][theta11][x[0],x[1],x[2],x[3]]
      + Derivative[1,0,1,0][theta12][x[0],x[1],x[2],x[3]]
      + Derivative[1,0,0,1][theta13][x[0],x[1],x[2],x[3]]
      - Derivative[1,1,0,0][theta11][x[0],x[1],x[2],x[3]]
      - Derivative[1,1,0,0][theta22][x[0],x[1],x[2],x[3]]
      - Derivative[1,1,0,0][theta33][x[0],x[1],x[2],x[3]]
   ))

  Print["component 01:"]
  Print[Series[einstein[0,1]-testeinstein[0,1],{eps,0,1}]]
  (* Abort[] *)

 testeinstein[1,2]:=(0 + eps(
 - 2 (2 Hp + H^2) Derivative[0,0,0,0][theta12][x[0],x[1],x[2],x[3]]
 - Derivative[0,1,1,0][h0][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,1,0][hD][x[0],x[1],x[2],x[3]]
 + Derivative[2,0,0,0][theta12][x[0],x[1],x[2],x[3]]
 + 2 H Derivative[1,0,0,0][theta12][x[0],x[1],x[2],x[3]]
 + 1/2 (
      + Derivative[1,0,1,0][h1][x[0],x[1],x[2],x[3]] 
      + Derivative[1,1,0,0][h2][x[0],x[1],x[2],x[3]] 
      )
 + H ( 
      + Derivative[0,0,1,0][h1][x[0],x[1],x[2],x[3]] 
      + Derivative[0,1,0,0][h2][x[0],x[1],x[2],x[3]] 
      )
 + Derivative[0,1,1,0][theta11][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,2,0][theta12][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,1,1][theta13][x[0],x[1],x[2],x[3]]
 + Derivative[0,2,0,0][theta12][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,1,0][theta22][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,0,1][theta23][x[0],x[1],x[2],x[3]]
 - (
        + Derivative[0,2,0,0][theta12][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta12][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta12][x[0],x[1],x[2],x[3]]
        )
 - (
        + Derivative[0,1,1,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,1,1,0][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[0,1,1,0][theta33][x[0],x[1],x[2],x[3]]
        )
 )
 ); 

  Print["component 12:"]
  Print[Series[einstein[1,2]-testeinstein[1,2],{eps,0,1}]]
 (* Abort[] *)

 testeinstein[1,1]:=( - (2 Hp + H^2) + eps(
 + 2 (2 Hp + H^2)(
      + Derivative[0,0,0,0][h0][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,0][hD][x[0],x[1],x[2],x[3]]
      - Derivative[0,0,0,0][theta11][x[0],x[1],x[2],x[3]]
      )
 + 2 H Derivative[1,0,0,0][h0][x[0],x[1],x[2],x[3]]
 + Derivative[0,2,0,0][h0][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,2,0][h0][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,0,2][h0][x[0],x[1],x[2],x[3]]
 - Derivative[0,2,0,0][h0][x[0],x[1],x[2],x[3]]
 - Derivative[0,2,0,0][hD][x[0],x[1],x[2],x[3]]
 - Derivative[0,0,2,0][hD][x[0],x[1],x[2],x[3]]
 - Derivative[0,0,0,2][hD][x[0],x[1],x[2],x[3]]
 + Derivative[0,2,0,0][hD][x[0],x[1],x[2],x[3]]
 + 2 Derivative[2,0,0,0][hD][x[0],x[1],x[2],x[3]]
 + 4 H Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
 + Derivative[2,0,0,0][theta11][x[0],x[1],x[2],x[3]]
 + 2 H Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
 + Derivative[1,1,0,0][h1][x[0],x[1],x[2],x[3]] 
 - Derivative[1,1,0,0][h1][x[0],x[1],x[2],x[3]] 
 - Derivative[1,0,1,0][h2][x[0],x[1],x[2],x[3]] 
 - Derivative[1,0,0,1][h3][x[0],x[1],x[2],x[3]] 
 + 2 H ( 
      + Derivative[0,1,0,0][h1][x[0],x[1],x[2],x[3]] 
      - Derivative[0,1,0,0][h1][x[0],x[1],x[2],x[3]] 
      - Derivative[0,0,1,0][h2][x[0],x[1],x[2],x[3]] 
      - Derivative[0,0,0,1][h3][x[0],x[1],x[2],x[3]] 
      )
 + 2(
   + Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]  
   + Derivative[0,1,1,0][theta12][x[0],x[1],x[2],x[3]]
   + Derivative[0,1,0,1][theta13][x[0],x[1],x[2],x[3]]
   )
 - (
        + Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta11][x[0],x[1],x[2],x[3]]
        )
 - (
      + Derivative[0,2,0,0][theta11][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,2,0][theta22][x[0],x[1],x[2],x[3]]
      + Derivative[0,0,0,2][theta33][x[0],x[1],x[2],x[3]]
      + 2 Derivative[0,1,1,0][theta12][x[0],x[1],x[2],x[3]]
      + 2 Derivative[0,1,0,1][theta13][x[0],x[1],x[2],x[3]]
      + 2 Derivative[0,0,1,1][theta23][x[0],x[1],x[2],x[3]]
        )
 + (
        - Derivative[2,0,0,0][theta11][x[0],x[1],x[2],x[3]]
        - Derivative[2,0,0,0][theta22][x[0],x[1],x[2],x[3]]
        - Derivative[2,0,0,0][theta33][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,2,0][theta33][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[0,0,0,2][theta33][x[0],x[1],x[2],x[3]]
        )
 - 2 H (
        + Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
        + Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
        + Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
       )
 )
 ); 

  Print["component 11:"]
  Print[Series[einstein[1,1]-testeinstein[1,1],{eps,0,1}]]
 (* Abort[] *)



