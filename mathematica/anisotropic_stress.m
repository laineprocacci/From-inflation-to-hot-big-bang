(*anisotropic stress*)

(*metric*)

 gdown[0,0]:=a[x[0]]^2(-1-2 eps h0[x[0],x[1],x[2],x[3]]);
 gdown[0,1]:=a[x[0]]^2(-eps h1[x[0],x[1],x[2],x[3]]);
 gdown[0,2]:=a[x[0]]^2(-eps h2[x[0],x[1],x[2],x[3]]);
 gdown[0,3]:=a[x[0]]^2(-eps h3[x[0],x[1],x[2],x[3]]);
 gdown[1,0]:=gdown[0,1];
 gdown[2,0]:=gdown[0,2];
 gdown[3,0]:=gdown[0,3];
 gdown[1,1]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]+2 eps theta11[x[0],x[1],x[2],x[3]]);
 gdown[2,2]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]+2 eps theta22[x[0],x[1],x[2],x[3]]);
 gdown[3,3]:=a[x[0]]^2(1-2 eps hD[x[0],x[1],x[2],x[3]]+2 eps theta33[x[0],x[1],x[2],x[3]]);
 gdown[1,2]:=a[x[0]]^2(2 eps theta12[x[0],x[1],x[2],x[3]]);
 gdown[1,3]:=a[x[0]]^2(2 eps theta13[x[0],x[1],x[2],x[3]]);
 gdown[2,3]:=a[x[0]]^2(2 eps theta23[x[0],x[1],x[2],x[3]]);
 gdown[2,1]:=gdown[1,2];
 gdown[3,1]:=gdown[1,3];
 gdown[3,2]:=gdown[2,3];

(*inverse metric*)

 gup[0,0]:=a[x[0]]^(-2)(-1+2 eps h0[x[0],x[1],x[2],x[3]]);
 gup[0,1]:=a[x[0]]^(-2)(-eps h1[x[0],x[1],x[2],x[3]]);
 gup[0,2]:=a[x[0]]^(-2)(-eps h2[x[0],x[1],x[2],x[3]]);
 gup[0,3]:=a[x[0]]^(-2)(-eps h3[x[0],x[1],x[2],x[3]]);
 gup[1,0]:=gup[0,1];
 gup[2,0]:=gup[0,2];
 gup[3,0]:=gup[0,3];
 gup[1,1]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]]-2 eps theta11[x[0],x[1],x[2],x[3]]);
 gup[2,2]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]]-2 eps theta22[x[0],x[1],x[2],x[3]]);
 gup[3,3]:=a[x[0]]^(-2)(1+2 eps hD[x[0],x[1],x[2],x[3]]-2 eps theta33[x[0],x[1],x[2],x[3]]);
 gup[1,2]:=a[x[0]]^(-2)(-2 eps theta12[x[0],x[1],x[2],x[3]]);
 gup[1,3]:=a[x[0]]^(-2)(-2 eps theta13[x[0],x[1],x[2],x[3]]);
 gup[2,3]:=a[x[0]]^(-2)(-2 eps theta23[x[0],x[1],x[2],x[3]]);
 gup[2,1]:=gup[1,2];
 gup[3,1]:=gup[1,3];
 gup[3,2]:=gup[2,3];

(*cristoffel symbols*)

 gamma[rho_,mu_,nu_]:=1/2 Sum[gup[rho,sigma](D[gdown[sigma,mu],x[nu]]+
  D[gdown[sigma,nu],x[mu]]-D[gdown[mu,nu],x[sigma]]),{sigma,0,3}]

(*background*)

  calH=a'[x[0]]/a[x[0]];
  calHp=a''[x[0]]/a[x[0]] - (a'[x[0]]/a[x[0]])^2;

(*4-velocity*)

  uup[0]:=a[x[0]]^(-1)( 1 - eps h0[x[0],x[1],x[2],x[3]]);
  uup[1]:=a[x[0]]^(-1)(  eps v1[x[0],x[1],x[2],x[3]]);
  uup[2]:=a[x[0]]^(-1)(  eps v2[x[0],x[1],x[2],x[3]]);
  uup[3]:=a[x[0]]^(-1)(  eps v3[x[0],x[1],x[2],x[3]]);

  udown[0]:=a[x[0]] ( -1 - eps h0[x[0],x[1],x[2],x[3]]);
  udown[1]:=a[x[0]] eps (v1[x[0],x[1],x[2],x[3]] - h1[x[0],x[1],x[2],x[3]]);
  udown[2]:=a[x[0]] eps (v2[x[0],x[1],x[2],x[3]] - h2[x[0],x[1],x[2],x[3]]);
  udown[3]:=a[x[0]] eps (v3[x[0],x[1],x[2],x[3]] - h3[x[0],x[1],x[2],x[3]]);

  testuu:=Sum[uup[alpha]udown[alpha],{alpha,0,3}];
  Print["test 4-velocity"]
  Print[Series[testuu,{eps,0,1}]]
  (* Abort[] *)

(*tensors of velocity derivatives*)

  udowndown[rho_,sigma_]:=D[udown[rho],x[sigma]]-
    Sum[gamma[alpha,rho,sigma] udown[alpha],{alpha,0,3}];
  uupdown[rho_,sigma_]:=D[uup[rho],x[sigma]]+
    Sum[gamma[rho,sigma,alpha] uup[alpha],{alpha,0,3}];
  utrace:=Sum[uupdown[alpha,alpha],{alpha,0,3}];

(*test tensors of velocity derivatives*)
 
  testudowndown[0,0]:=0;
  testudowndown[0,1]:=-eps a[x[0]] ( 
                          + calH Derivative[0,0,0,0][v1][x[0],x[1],x[2],x[3]]
			  );
  testudowndown[1,0]:= eps a[x[0]] (   
                          + Derivative[1,0,0,0][v1][x[0],x[1],x[2],x[3]]
                          - Derivative[1,0,0,0][h1][x[0],x[1],x[2],x[3]]
                          - calH Derivative[0,0,0,0][h1][x[0],x[1],x[2],x[3]]
                          + Derivative[0,1,0,0][h0][x[0],x[1],x[2],x[3]]
			     );
  testudowndown[1,1]:=a[x[0]] calH + 1/2 eps a[x[0]] (   
                   - 2 calH Derivative[0,0,0,0][h0][x[0],x[1],x[2],x[3]]
                        + 2 Derivative[0,1,0,0][v1][x[0],x[1],x[2],x[3]]
                   - 4 calH Derivative[0,0,0,0][hD][x[0],x[1],x[2],x[3]]
                        - 2 Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
                   + 4 calH Derivative[0,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                        + 2 Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
			     );
  testudowndown[1,2]:=1/2 eps a[x[0]] (   
                          + 2 Derivative[0,0,1,0][v1][x[0],x[1],x[2],x[3]]
                          - Derivative[0,0,1,0][h1][x[0],x[1],x[2],x[3]]
                          + Derivative[0,1,0,0][h2][x[0],x[1],x[2],x[3]]
                   + 4 calH Derivative[0,0,0,0][theta12][x[0],x[1],x[2],x[3]]
                        + 2 Derivative[1,0,0,0][theta12][x[0],x[1],x[2],x[3]]
			     );
  testudowndown[2,1]:=1/2 eps a[x[0]] (   
                          + 2 Derivative[0,1,0,0][v2][x[0],x[1],x[2],x[3]]
                          + Derivative[0,0,1,0][h1][x[0],x[1],x[2],x[3]]
                          - Derivative[0,1,0,0][h2][x[0],x[1],x[2],x[3]]
                   + 4 calH Derivative[0,0,0,0][theta12][x[0],x[1],x[2],x[3]]
                        + 2 Derivative[1,0,0,0][theta12][x[0],x[1],x[2],x[3]]
			     );

  testutrace:=3 a[x[0]]^(-1) calH + a[x[0]]^(-1) eps (
                   - 3 calH Derivative[0,0,0,0][h0][x[0],x[1],x[2],x[3]]
                         +  Derivative[0,1,0,0][v1][x[0],x[1],x[2],x[3]]
                         +  Derivative[0,0,1,0][v2][x[0],x[1],x[2],x[3]]
                         +  Derivative[0,0,0,1][v3][x[0],x[1],x[2],x[3]]
                        - 3 Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
                         +  Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                         +  Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
                         +  Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
                             );

  Print["test velocity derivatives"]
  Print[Series[udowndown[0,0]-testudowndown[0,0],{eps,0,1}]]
  Print[Series[udowndown[0,1]-testudowndown[0,1],{eps,0,1}]]
  Print[Series[udowndown[1,0]-testudowndown[1,0],{eps,0,1}]]
  Print[Series[udowndown[1,1]-testudowndown[1,1],{eps,0,1}]]
  Print[Series[udowndown[1,2]-testudowndown[1,2],{eps,0,1}]]
  Print[Series[udowndown[2,1]-testudowndown[2,1],{eps,0,1}]]

  Print["test velocity divergence"]
  Print[Series[utrace-testutrace,{eps,0,1}]]
  (* Abort[] *)

(*transverse velocity projectors*)

  projdowndown[mu_,nu_]:=gdown[mu,nu]+udown[mu] udown[nu];
  projdownup[mu_,nu_]:=Sum[gdown[mu,alpha1] gup[alpha1,nu],{alpha1,0,3}]+
                        udown[mu] uup[nu]; 


  testprojdowndown[0,0]:=0;
  testprojdowndown[0,1]:=a[x[0]]^2 eps(
                     - Derivative[0,0,0,0][v1][x[0],x[1],x[2],x[3]]
                          );
  testprojdowndown[1,0]:=testprojdowndown[0,1];
  testprojdowndown[1,2]:=a[x[0]]^2 eps (
                     + 2 Derivative[0,0,0,0][theta12][x[0],x[1],x[2],x[3]]
			 );     
  testprojdowndown[2,1]:=testprojdowndown[1,2];
testprojdowndown[1,1]:=a[x[0]]^2 (1 + eps (
                     - 2     Derivative[0,0,0,0][hD][x[0],x[1],x[2],x[3]]
                     + 2     Derivative[0,0,0,0][theta11][x[0],x[1],x[2],x[3]]
		     ));  

  Print["test velocity transverse projector"]
  Print[Series[projdowndown[0,0]-testprojdowndown[0,0],{eps,0,1}]]
  Print[Series[projdowndown[0,1]-testprojdowndown[0,1],{eps,0,1}]]
  Print[Series[projdowndown[1,0]-testprojdowndown[1,0],{eps,0,1}]]
  Print[Series[projdowndown[1,2]-testprojdowndown[1,2],{eps,0,1}]]
  Print[Series[projdowndown[2,1]-testprojdowndown[2,1],{eps,0,1}]]
  Print[Series[projdowndown[1,1]-testprojdowndown[1,1],{eps,0,1}]]
  (* Abort[] *)

(*shear part*)

  shear[mu_,nu_]:=-Sum[projdownup[mu,rho1] projdownup[nu,sigma1](
		  +udowndown[rho1,sigma1]
                  +udowndown[sigma1,rho1]
                  -2/3 gdown[rho1,sigma1] utrace)
                   ,{rho1,0,3},{sigma1,0,3}]

  testshear[0,0]:=0;
  testshear[0,1]:=0;
  testshear[1,0]:=0;
  testshear[1,2]:=-a[x[0]] eps (
                         +  Derivative[0,0,1,0][v1][x[0],x[1],x[2],x[3]]
                         +  Derivative[0,1,0,0][v2][x[0],x[1],x[2],x[3]]
                        + 2 Derivative[1,0,0,0][theta12][x[0],x[1],x[2],x[3]]
			 );     
  testshear[2,1]:=testshear[1,2];
  testshear[1,1]:=-a[x[0]] eps (
                        + 2 Derivative[0,1,0,0][v1][x[0],x[1],x[2],x[3]]
			- 2/3 (  Derivative[0,1,0,0][v1][x[0],x[1],x[2],x[3]]
                              +  Derivative[0,0,1,0][v2][x[0],x[1],x[2],x[3]]
                              +  Derivative[0,0,0,1][v3][x[0],x[1],x[2],x[3]]
				)
                        + 2 Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
			- 2/3 (  Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                              +  Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
                              +  Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
				)
			 );     

  Print["test shear part"]
  Print[Series[shear[0,0]-testshear[0,0],{eps,0,1}]]
  Print[Series[shear[0,1]-testshear[0,1],{eps,0,1}]]
  Print[Series[shear[1,0]-testshear[1,0],{eps,0,1}]]
  Print[Series[shear[1,2]-testshear[1,2],{eps,0,1}]]
  Print[Series[shear[2,1]-testshear[2,1],{eps,0,1}]]
  Print[Series[shear[1,1]-testshear[1,1],{eps,0,1}]]
  (* Abort[] *)

(*bulk part*)

  bulk[mu_,nu_]:=- ( projdowndown[mu,nu] utrace)

  testbulk[0,0]:=0;
  testbulk[0,1]:=-3 a[x[0]] eps(
                     -calH Derivative[0,0,0,0][v1][x[0],x[1],x[2],x[3]]
                          );
  testbulk[1,0]:=testbulk[0,1];
  testbulk[1,2]:=-3 a[x[0]] eps (
                     + 2 calH Derivative[0,0,0,0][theta12][x[0],x[1],x[2],x[3]]
			 );     
  testbulk[2,1]:=testbulk[1,2];
  testbulk[1,1]:=-3 a[x[0]] calH - 3 a[x[0]] eps (
                     + 2 calH Derivative[0,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                     -   calH Derivative[0,0,0,0][h0][x[0],x[1],x[2],x[3]]
                     -        Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
                     - 2 calH Derivative[0,0,0,0][hD][x[0],x[1],x[2],x[3]]
                     +1/3 (  Derivative[0,1,0,0][v1][x[0],x[1],x[2],x[3]]
                         +   Derivative[0,0,1,0][v2][x[0],x[1],x[2],x[3]]
                         +   Derivative[0,0,0,1][v3][x[0],x[1],x[2],x[3]]
			     )
                     +1/3 (  Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                         +   Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
                         +   Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
			     )
			 );     

  Print["test bulk part"]
  Print[Series[bulk[0,0]-testbulk[0,0],{eps,0,1}]]
  Print[Series[bulk[0,1]-testbulk[0,1],{eps,0,1}]]
  Print[Series[bulk[1,0]-testbulk[1,0],{eps,0,1}]]
  Print[Series[bulk[1,2]-testbulk[1,2],{eps,0,1}]]
  Print[Series[bulk[2,1]-testbulk[2,1],{eps,0,1}]]
  Print[Series[bulk[1,1]-testbulk[1,1],{eps,0,1}]]
  (* Abort[] *)
  
