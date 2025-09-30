
(*bose distribution*)

 nb[x_]:=1/(Exp[x]-1)

(*photon number density per degree of freedom*)

 ngamma:=Integrate[p^2/(2Pi^2) nb[p],{p,0,Infinity}]

 ngammatest:=Zeta[3]/Pi^2;
 
 Print[FullSimplify[ngamma-ngammatest]]  
 Print["number density checked"]

(*photon energy density per degree of freedom*)

 egamma:=Integrate[p^3/(2Pi^2) nb[p],{p,0,Infinity}]

 egammatest:=Pi^2/30;
 
 Print[FullSimplify[egamma-egammatest]]  
 Print["energy density checked"]





