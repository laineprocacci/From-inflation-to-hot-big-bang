
(*energy-momentum conservation for varphi + fluid up to first order in perturbations*)

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

  calH=a'[x[0]]/a[x[0]];
  calHp=a''[x[0]]/a[x[0]] - (a'[x[0]]/a[x[0]])^2;

(*cristoffel symbols*)

 gamma[rho_,mu_,nu_]:=1/2 Sum[gup[rho,sigma](D[gdown[sigma,mu],x[nu]]+
  D[gdown[sigma,nu],x[mu]]-D[gdown[mu,nu],x[sigma]]),{sigma,0,3}]

(*energy-momentum tensor for varphi plus perfect fluid*)

 t[0,0]:=(barphi'[x[0]])^2/2+
    eps barphi'[x[0]] Derivative[1,0,0,0][dphi][x[0],x[1],x[2],x[3]]+
    a[x[0]]^2 (bare[x[0]]+eps(deltae[x[0],x[1],x[2],x[3]]
   +2 bare[x[0]] h0[x[0],x[1],x[2],x[3]]));
 t[0,1]:=eps(barphi'[x[0]] Derivative[0,1,0,0][dphi][x[0],x[1],x[2],x[3]]
  -(barphi'[x[0]])^2/2 h1[x[0],x[1],x[2],x[3]] )+
   a[x[0]]^2 eps(-bare[x[0]](v1[x[0],x[1],x[2],x[3]]-h1[x[0],x[1],x[2],x[3]])
  -barp[x[0]] v1[x[0],x[1],x[2],x[3]]);
 t[0,2]:=eps(barphi'[x[0]] Derivative[0,0,1,0][dphi][x[0],x[1],x[2],x[3]]
  -(barphi'[x[0]])^2/2 h2[x[0],x[1],x[2],x[3]] )+
   a[x[0]]^2 eps(-bare[x[0]](v2[x[0],x[1],x[2],x[3]]-h2[x[0],x[1],x[2],x[3]])
  -barp[x[0]] v2[x[0],x[1],x[2],x[3]]);
 t[0,3]:=eps(barphi'[x[0]] Derivative[0,0,0,1][dphi][x[0],x[1],x[2],x[3]]
  -(barphi'[x[0]])^2/2 h3[x[0],x[1],x[2],x[3]] )+
   a[x[0]]^2 eps(-bare[x[0]](v3[x[0],x[1],x[2],x[3]]-h3[x[0],x[1],x[2],x[3]])
  -barp[x[0]] v3[x[0],x[1],x[2],x[3]]);
 t[1,0]:=t[0,1];
 t[2,0]:=t[0,2];
 t[3,0]:=t[0,3];
 t[1,1]:=(barphi'[x[0]])^2/2+
  eps(2(hD[x[0],x[1],x[2],x[3]]-theta11[x[0],x[1],x[2],x[3]])(-(barphi'[x[0]])^2/2)
  -h0[x[0],x[1],x[2],x[3]] (barphi'[x[0]])^2
  +barphi'[x[0]] Derivative[1,0,0,0][dphi][x[0],x[1],x[2],x[3]])+
   a[x[0]]^2 (barp[x[0]]+eps(deltap[x[0],x[1],x[2],x[3]]
  -2 barp[x[0]] hD[x[0],x[1],x[2],x[3]]
  +2 barp[x[0]] theta11[x[0],x[1],x[2],x[3]]
  +             pi11[x[0],x[1],x[2],x[3]]));
 t[2,2]:=(barphi'[x[0]])^2/2+
  eps(2(hD[x[0],x[1],x[2],x[3]]-theta22[x[0],x[1],x[2],x[3]])(-(barphi'[x[0]])^2/2)
  -h0[x[0],x[1],x[2],x[3]] (barphi'[x[0]])^2
  +barphi'[x[0]] Derivative[1,0,0,0][dphi][x[0],x[1],x[2],x[3]])+
  a[x[0]]^2(barp[x[0]]+eps(deltap[x[0],x[1],x[2],x[3]]
  -2 barp[x[0]] hD[x[0],x[1],x[2],x[3]]
  +2 barp[x[0]] theta22[x[0],x[1],x[2],x[3]]
  +             pi22[x[0],x[1],x[2],x[3]]));
 t[3,3]:=(barphi'[x[0]])^2/2+  
  eps(2(hD[x[0],x[1],x[2],x[3]]-theta33[x[0],x[1],x[2],x[3]])(-(barphi'[x[0]])^2/2)
  -h0[x[0],x[1],x[2],x[3]] (barphi'[x[0]])^2
  +barphi'[x[0]] Derivative[1,0,0,0][dphi][x[0],x[1],x[2],x[3]])+
   a[x[0]]^2(barp[x[0]]+eps(deltap[x[0],x[1],x[2],x[3]]
  -2 barp[x[0]] hD[x[0],x[1],x[2],x[3]]
  +2 barp[x[0]] theta33[x[0],x[1],x[2],x[3]]
  +             pi33[x[0],x[1],x[2],x[3]]));
 t[1,2]:=eps theta12[x[0],x[1],x[2],x[3]] (barphi'[x[0]])^2+ 
  a[x[0]]^2 eps(2 barp[x[0]] theta12[x[0],x[1],x[2],x[3]]+pi12[x[0],x[1],x[2],x[3]]);
 t[1,3]:=eps theta13[x[0],x[1],x[2],x[3]] (barphi'[x[0]])^2+ 
  a[x[0]]^2 eps(2 barp[x[0]] theta13[x[0],x[1],x[2],x[3]]+pi13[x[0],x[1],x[2],x[3]]);
 t[2,3]:=eps theta23[x[0],x[1],x[2],x[3]] (barphi'[x[0]])^2+ 
  a[x[0]]^2 eps(2 barp[x[0]] theta23[x[0],x[1],x[2],x[3]]+pi23[x[0],x[1],x[2],x[3]]);
 t[2,1]:=t[1,2]; 
 t[3,1]:=t[1,3];
 t[3,2]:=t[2,3];

(* Print["00: ",t[0,0]];
 Print["01: ",t[0,1]];
 Print["11: ",t[1,1]];
 Print["12: ",t[1,2]]; *)

(*energy-momentum conservation equation*)

 conserve[nu_]:=Sum[D[t[mu,nu],x[alpha]]gup[alpha,mu],{alpha,0,3},{mu,0,3}]-
  Sum[(t[mu,beta]gamma[beta,nu,alpha]+
  t[beta,nu]gamma[beta,mu,alpha])gup[alpha,mu],{alpha,0,3},{mu,0,3},{beta,0,3}]

(*test energy-momentum conservation*)

(*component nu=0*)

  testconserve[0]:=-(bare'[x[0]]+3 calH(bare[x[0]]+barp[x[0]]))-
   1/a[x[0]]^2 barphi'[x[0]](barphi''[x[0]]+2calH barphi'[x[0]])+
   eps(
   +1/a[x[0]]^2(
   + barphi'[x[0]] (-Derivative[2,0,0,0][dphi][x[0],x[1],x[2],x[3]]
                    +Derivative[0,2,0,0][dphi][x[0],x[1],x[2],x[3]]
                    +Derivative[0,0,2,0][dphi][x[0],x[1],x[2],x[3]]
                    +Derivative[0,0,0,2][dphi][x[0],x[1],x[2],x[3]]
   )
   - ( barphi''[x[0]] 
   + 4 calH barphi'[x[0]] ) Derivative[1,0,0,0][dphi][x[0],x[1],x[2],x[3]]
   + barphi'[x[0]]^2 (+Derivative[1,0,0,0][h0][x[0],x[1],x[2],x[3]]
                      +3 Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
                      -Derivative[0,1,0,0][h1][x[0],x[1],x[2],x[3]]
                      -Derivative[0,0,1,0][h2][x[0],x[1],x[2],x[3]]
                      -Derivative[0,0,0,1][h3][x[0],x[1],x[2],x[3]]
                      -Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                      -Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
                      -Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
   )
   +2 barphi'[x[0]] ( 
   + barphi''[x[0]] + 2 calH barphi'[x[0]] ) h0[x[0],x[1],x[2],x[3]]
   )
  -Derivative[1,0,0,0][deltae][x[0],x[1],x[2],x[3]]
  -3 calH (+ Derivative[0,0,0,0][deltae][x[0],x[1],x[2],x[3]]
           + Derivative[0,0,0,0][deltap][x[0],x[1],x[2],x[3]] 
          )
  +3 (bare[x[0]]+barp[x[0]]) Derivative[1,0,0,0][hD][x[0],x[1],x[2],x[3]]
  - (bare[x[0]]+barp[x[0]]) ( Derivative[0,1,0,0][v1][x[0],x[1],x[2],x[3]]
                             +Derivative[0,0,1,0][v2][x[0],x[1],x[2],x[3]]
                             +Derivative[0,0,0,1][v3][x[0],x[1],x[2],x[3]]
                             +Derivative[1,0,0,0][theta11][x[0],x[1],x[2],x[3]]
                             +Derivative[1,0,0,0][theta22][x[0],x[1],x[2],x[3]]
                             +Derivative[1,0,0,0][theta33][x[0],x[1],x[2],x[3]]
                            )
  -calH ( Derivative[0,0,0,0][pi11][x[0],x[1],x[2],x[3]]
         +Derivative[0,0,0,0][pi22][x[0],x[1],x[2],x[3]]
         +Derivative[0,0,0,0][pi33][x[0],x[1],x[2],x[3]])
  );

 Print[Series[conserve[0]-testconserve[0],{eps,0,1}]]
 Print["energy-momentum conservation tested: nu=0"]

(*component nu=1*)

 testconserve[1]:=eps(
 - 1/a[x[0]]^2 (barphi''[x[0]] + 2 calH barphi'[x[0]]) Derivative[0,1,0,0][dphi][x[0],x[1],x[2],x[3]]
 + Derivative[0,1,0,0][deltap][x[0],x[1],x[2],x[3]]
 + (bare[x[0]]+barp[x[0]]) (  Derivative[0,1,0,0][h0][x[0],x[1],x[2],x[3]]
                            + Derivative[1,0,0,0][v1][x[0],x[1],x[2],x[3]]
                            - Derivative[1,0,0,0][h1][x[0],x[1],x[2],x[3]]
			    + 4 calH ( Derivative[0,0,0,0][v1][x[0],x[1],x[2],x[3]]
                                     - Derivative[0,0,0,0][h1][x[0],x[1],x[2],x[3]]
				       )
                             )
 + (bare'[x[0]]+barp'[x[0]]) (  Derivative[0,0,0,0][v1][x[0],x[1],x[2],x[3]]
                            -   Derivative[0,0,0,0][h1][x[0],x[1],x[2],x[3]]
                            )
 + Derivative[0,1,0,0][pi11][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,1,0][pi12][x[0],x[1],x[2],x[3]]
 + Derivative[0,0,0,1][pi13][x[0],x[1],x[2],x[3]]
  );

 Print[Series[conserve[1]-testconserve[1],{eps,0,1}]]
 Print["energy-momentum conservation tested: nu=1"]


