
(* --- equation for dot{S_T} from sec.7.5 --- *)

Print[" equation for dot{S_T} from sec.7.5 "]

(* 7.83 / nu=0 of energy-momentum conservation *)

  zero[1] = 1/a^2 (
		   (2 a Ups barphip + a^2 Vphi ) dphip
		 +  a barphip^2 dUps
		 +  a^2 barphip (dVphi - rho )
	         -  a Ups barphip^2 h0 ) + (
                 -  dep
                 -  3 H (de + dp)  
                 +  (bare + barp) (comb + nabla2 vmh - h0p )
		 )
		     
(* insert perturbations in terms of gauge-invariant quantities / without Y for now *)

   subst1   = { dphip -> 1/a(Qp - H Q) (* 7.62 *),
                dUps -> Upsp/(a barphip) Q - a UpsT/(H bareT) ST (* 7.64 *),
                dVphi -> Vphip/(a barphip) Q - a VphiT/(H bareT) ST (* 7.65 *),
                h0 -> (H-Hp/H)/(a barphip) Q - 4 Pi G a^2/H^2 Sv (* 7.68 *),
                dep -> barep/(a barphip) Qp + ( barepp/(a barphip) -
		        barep/(a barphip)(H + barphipp/barphip)) Q -
                        a/H STp-a/H(H-Hp/H) ST (* 7.90 *),
                de -> barep/(a barphip) Q -a/H ST (* 7.66 *),
                dp -> barpp/(a barphip) Q -a barpT/(H bareT) ST (* 7.66 *),
                comb -> -2 (2H+Hp/H) h0 + 4 Pi G a^2/H (dp - de + 2/3 nabla2 Pi) (* 7.55 *),
                vmh -> Q/(a barphip) - Sv/H/(bare+barp) (* 7.67 *),
                h0p -> (H - Hp/H)/(a barphip) Qp +(
			2 Hp - H^2 - Hpp/H + Hp^2/H^2 - (H-Hp/H)barphipp/barphip
                        )/(a barphip) Q - 4 Pi G a^2/H^2 ( Svp + 2 (H - Hp/H) Sv ) (* 7.92 *),
                barepp -> -3 Hp (bare+barp) - 3 H (barep + barpp) + 
		           (Upsp - H Ups) barphip^2/a + 2 Ups barphip barphipp/a +
                           barphipp Vphi + barphip Vphip (* 7.94 *),
                Hpp -> -2 H Hp + 4 H^3 + 4 Pi G a^2 (barep - barpp) (* 7.73 *),
                Svp -> a barpT / bareT ST - (3 H + 4 Pi G barphip^2 / H) Sv +
		(bare+barp)H/(a barphip) Qp - (bare+barp)H/(a barphip) (
		barphipp/barphip + H - Hp/H) Q -2/3 H nabla2 Pi (* 7.87 *),
                Qp -> Qpcomb + (barphipp/barphip + H - Hp/H) Q,
                barep -> -3 H(bare+barp) - (barphipp + 2 H barphip)barphip/a^2 (* 1.72 *),
                barphipp -> -(2 H + a Ups) barphip - a^2 Vphi (* 1.66 *)
     }

     zero[2]=Collect[zero[1]//.subst1,{Qpcomb,Q,ST,Sv,Pi,rho},FullSimplify]

(* insert 7.95 in conformal time *)

     testzero = rho barphip H / a^2 + 
		( Ups barphip^2/a^2 + 8 Pi G a bare (bare + barp)/H )(Rp / a - 4 Pi G a/ H Sv) + 
		nabla2/a^2 ((bare + barp) R + Sv) -
                STp/a - 
		(3 H/a(1+barpT/bareT)+4 Pi a G (bare+barp)/H - (Hp - H^2)/H/a - 
		 barphip (UpsT barphip + a VphiT)/bareT/a^2) ST

     subst2 = {R -> -H/(a barphip) Q,Rp->-H/(a barphip) Qpcomb} (* 7.76 *)

     zero[3]=Collect[(-H/a^2)zero[2]-(testzero/.subst2),{Qpcomb,Q,ST,Sv,Pi,rho},FullSimplify]

(* identify remainder *)

     remainder = Qpcomb/(a^3 barphip)(bare+barp)(-2 H^2 - Hp + 4 Pi G a^2(bare-barp))+
		 Sv 4 Pi G (bare + barp)/H^2 (-3 H^2 + 4 Pi G (2 a^2 bare + barphip^2))

(* check that all terms get subtracted *)

     zero[4]=Collect[zero[3]-remainder,{Qpcomb,Sv},FullSimplify]
		     
     Print[zero[4]]

(**************************************************)
		     
(* --- equation for dot{R_varphi} from sec.8.4 --- *)

  Print[" equation for dot{R_varphi} from sec.8.4 "]

(* 8.23 *)

  zero[1] =     (
                   dphihatpp-nabla2 dphihat
		 + (a^2 Vphiphi - Hp - H^2) dphihat
                 + a^2 Vphichi dchihat
                 - comb a barphip + 2 h0 a^3 Vphi 
		 )
		     
(* insert perturbations in terms of gauge-invariant quantities *)

   subst1   = { dphihat -> Qphi - preFphi Y (* 8.30 *), 
                dphihatp -> Qphip - preFphip Y - preFphi Yp, 
                dphihatpp -> Qphipp - preFphipp Y - 2 preFphip Yp - preFphi Ypp, 
                dchihat -> Qchi - preFchi Y (* 8.30 *), 
                preFphi -> a barphip/H, 
                preFphip -> a barphip (-1-Hp/H^2) - a^3 Vphi/H (* 8.30 *), 
                preFphipp -> a barphip/H ( H^2 + Hp - (2 + Hp/H^2) 8 Pi G (barphip^2+barchip^2)
			  -16 Pi G a^2 Vphi barphip/H - a^2 Vphiphi ) - a barchip/H (
                          8 Pi G a^2 (Vphi barchip + Vchi barphip)/H + a^2 Vphichi ) (* 8.31 *),
                preFchi -> a barchip/H, 
                h0 -> - Yp/H + 4 Pi G/(a H)( barphip dphihat + barchip dchihat ) (* 8.27 *),
                comb -> -2 (2H+Hp/H) h0 -(Ypp + 2 H Yp - nabla2 Y)/H - 8 Pi G a/H(
                         Vphi dphihat + Vchi dchihat ) (* 8.28 *)
     }

     zero[2]=Collect[zero[1]//.subst1,{Y,Yp,Ypp,Qphi,Qphip,Qphipp},FullSimplify]

(* insert 8.33 in conformal time *)

     testzero = Qphipp - nabla2 Qphi + mmphi Qphi + ( a^2 Vphichi
	        + (2 + Hp/H^2)8 Pi G barphip barchip
                + 8 Pi G a^2 (Vphi barchip + Vchi barphip)/H
		)(Qchi - barchip/barphip Qphi)

     subst2 = {mmphi->-H/(a barphip) preFphipp} (* 8.32 *)

(* check that all terms get subtracted *)
		     
     zero[3]=Collect[(zero[2]-(testzero/.subst2))/.subst1,
	 {Qphi,Qphip,Qphipp,Qchi,Qchip,Qchipp},FullSimplify]

     Print[zero[3]]

(**************************************************)

(* --- equation for dot{R_T} from sec.9.2 --- *)

Print[" equation for dot{R_T} from sec.9.2 "]

(* 7.83 / nu=0 of energy-momentum conservation *)

  zero[1] =  (
                 -  dep
                 -  3 H (de + dp)  
                 +  (bare + barp) (comb + nabla2 vmh - h0p )
		 )
		     
(* insert perturbations in terms of gauge-invariant quantities and Y *)

   subst1   = { h0 -> -Yp/H + 4 Pi G a^2/H (bare+barp) vmh (* 9.3 *),
                vmh -> - hmv,
                hmv -> Rv/H + Y/H (* 9.4 *),
                de -> - barep(RT/H+Y/H) (* 9.5 *),
                dp -> - barpp(RT/H+Y/H) (* 9.5 *),
                dep -> - barepp(RT/H+Y/H) - barep(RTp/H+Yp/H) + barep Hp (RT+Y)/H^2 (* 9.5 *),
                comb -> -2 (2H+Hp/H) h0 - (Ypp + 2 H Yp - nabla2 Y) /H + 
                        4 Pi G a^2/H (dp - de + 2/3 nabla2 Pi) (* 9.10 *),
                h0p -> (Hp Yp - H Ypp)/H^2 + (1-Hp/H^2)(-Rvp-Yp) +(
			2 Hp^2/H^3 - Hpp/H^2) H vmh  (* 9.11 *),
                Hp -> H^2 - 4 Pi G a^2 (bare + barp) (* 9.7 *),
                Hpp -> -2 H Hp + 4 H^3 + 4 Pi G a^2 (barep - barpp) (* 7.75 *),
                barep -> -3 H(bare+barp) (* 9.7 *),
                barepp -> -3 Hp (bare+barp) - 3 H (barep + barpp) (* 9.7 *),
                Rvp -> -2/3 H nabla2 Pi/(bare + barp) - barpp/(bare+barp)(Rv-RT) (* 9.8 *)
     }

     zero[2]=Collect[zero[1]//.subst1,{Y,Yp,Ypp,Rv,Rvp,RT,RTp},FullSimplify]

(* insert 7.95 in conformal time *)

     testzero = (bare + barp)(-3 RTp - nabla2 Rv/H + 4 Pi G a^2 barep/H^2 (RT-Rv))

     zero[3]=Collect[ (zero[2]-testzero)/.subst1,{Rv,Rvp,RT,RTp},FullSimplify]

     Print[zero[3]]

(**************************************************)

(* --- equation for dot{R_v} from sec.9.5 --- *)

Print[" equation for dot{R_v} from sec.9.5 "]

(* 9.49 / nu=i of energy-momentum conservation *)

  zero[1] =  (
                 +  dp1
                 +  (bare1 + barp1) h0  
                 +  (bare1p + barp1p) hmv1 + (bare1 + barp1) ( hmv1p + 4 H hmv1) 
                 +  2/3 nabla2 Pi1
		 ) (* 9.40 *)
		     
(* insert perturbations in terms of gauge-invariant quantities and Y *)

   subst1   = { h0 -> -Yp/H +(Hp/H^2 -1)Y - 4 Pi G a^2/H^2 (
		   (bare1+barp1) Rv1 + (bare2+barp2) Rv2 ) (* 9.45 *),
                vmh1 -> - hmv1,
                hmv1 -> (Rv1 + Y)/H (* 9.41 *),
                hmv1p -> (Rv1p + Yp)/H - (Rv1 + Y) Hp/H^2 (* 9.41 *),
                de1 -> - bare1p(RT1/H+Y/H) (* 9.42 *),
                dp1 -> - barp1p(RT1/H+Y/H) (* 9.41 *),
                de1p -> 3 (bare1p + barp1p) (RT1 + Y) + 3 (bare1+barp1) (RT1p + Yp)  (* 9.42 *),
                comb -> -2 (2H+Hp/H) h0 - (Ypp + 2 H Yp - nabla2 Y) /H + 
  		        4 Pi G a^2/H (dp1 + dp2 - de1 -de2 + 2/3 nabla2 (Pi1 + Pi2)) (* 9.46 *),
                h0p -> - Ypp/H^2  + (2Hp /H^2 - 1) Yp + (Hpp/H^2 - 2 Hp^2/H^3) Y +2 (
			 Hp/H - H) 4 Pi G a^2/H^2 ( (bare1+barp1) Rv1 + (bare2+barp2) Rv2 ) - 
		         4 Pi G a^2/H^2( (bare1+barp1) Rv1p + (bare2+barp2) Rv2p + 
                                (bare1p+barp1p) Rv1 + (bare2p+barp2p) Rv2 ) (* 9.47 *),
                Hp -> H^2 - 4 Pi G a^2 (bare1 + barp1 + bare2 + barp2) (* 9.43 *),
                Hpp -> -2 H Hp + 4 H^3 + 4 Pi G a^2 (bare1p + bare2p - barp1p - barp2p) (* 9.44 *),
                bare1p -> -3 H(bare1+barp1) (* 9.38 *),
                bare2p -> -3 H(bare2+barp2) (* 9.38 *),
                bare1pp -> -3 Hp (bare1+barp1) - 3 H (bare1p + barp1p) (* 9.38 *),
                bare2pp -> -3 Hp (bare2+barp2) - 3 H (bare2p + barp2p) (* 9.38 *)
                (* Rvp -> -2/3 H nabla2 Pi/(bare + barp) - barpp/(bare+barp)(Rv-RT) 9.8 *)
     }

     zero[2]=Collect[zero[1]//.subst1,{Y,Yp,Ypp,Rv1,Rv1p,RT1,RT1p,Rv2,Rv2p,RT2,RT2p},FullSimplify]

(* insert 9.49 *)

     testzero = 1/H(barp1p(Rv1-RT1)-4 Pi G a^2/H(bare1+barp1)(bare2+barp2)(Rv2-Rv1)+
                    (bare1+barp1)Rv1p + 2H/3 nabla2 Pi1)

     zero[3]=Collect[ (zero[2]-testzero)/.subst1,{Rv1,Rv1p,Rv2,Rv2p,RT1,RT1p},FullSimplify]

     Print[zero[3]]

(**************************************************)

(* --- equation for dot{R_T} from sec.9.5 --- *)

Print[" equation for dot{R_T} from sec.9.5 "]

(* 9.51 / nu=0 of energy-momentum conservation *)

  zero[1] =  (
                 -  de1p
                 -  3 H (de1 + dp1)   
                 +  (bare1 + barp1) ( comb + nabla2 vmh1 - h0p) 
		 )   (* 9.39 *)
		     
(* insert perturbations in terms of gauge-invariant quantities and Y *)

   subst1   = { h0 -> -Yp/H +(Hp/H^2 -1)Y - 4 Pi G a^2/H^2 (
		   (bare1+barp1) Rv1 + (bare2+barp2) Rv2 ) (* 9.45 *),
                vmh1 -> - hmv1,
                hmv1 -> (Rv1 + Y)/H (* 9.41 *),
                hmv1p -> (Rv1p + Yp)/H - (Rv1 + Y) Hp/H^2 (* 9.41 *),
                dp1 -> - barp1p(RT1+Y)/H (* 9.41 *),
                dp2 -> - barp2p(RT2+Y)/H (* 9.41 *),
                de1 -> 3 (bare1 + barp1)(RT1+Y) (* 9.42 *),
                de2 -> 3 (bare2 + barp2)(RT2+Y) (* 9.42 *),
                de1p -> 3 (bare1p + barp1p) (RT1 + Y) + 3 (bare1+barp1) (RT1p + Yp)  (* 9.42 *),
                comb -> -2 (2H+Hp/H) h0 - (Ypp + 2 H Yp - nabla2 Y) /H + 
  		        4 Pi G a^2/H (dp1 + dp2 - de1 -de2 + 2/3 nabla2 (Pi1 + Pi2)) (* 9.46 *),
                h0p -> - Ypp/H  + (2Hp /H^2 - 1) Yp + (Hpp/H^2 - 2 Hp^2/H^3) Y +2 (
			 Hp/H - H) 4 Pi G a^2/H^2 ( (bare1+barp1) Rv1 + (bare2+barp2) Rv2 ) - 
		         4 Pi G a^2/H^2( (bare1+barp1) Rv1p + (bare2+barp2) Rv2p + 
                                (bare1p+barp1p) Rv1 + (bare2p+barp2p) Rv2 ) (* 9.47 *),
                Hp -> H^2 - 4 Pi G a^2 (bare1 + barp1 + bare2 + barp2) (* 9.43 *),
                Hpp -> -2 H Hp + 4 H^3 + 4 Pi G a^2 (bare1p + bare2p - barp1p - barp2p) (* 9.44 *),
                bare1p -> -3 H(bare1+barp1) (* 9.38 *),
                bare2p -> -3 H(bare2+barp2) (* 9.38 *),
                bare1pp -> -3 Hp (bare1+barp1) - 3 H (bare1p + barp1p) (* 9.38 *),
                bare2pp -> -3 Hp (bare2+barp2) - 3 H (bare2p + barp2p) (* 9.38 *), 
                Rv1p -> -1/(bare1+barp1)(
                    barp1p(Rv1-RT1)-4 Pi G a^2/H(bare1+barp1)(bare2+barp2)(Rv2-Rv1)+
                    + 2H/3 nabla2 Pi1 ),
                Rv2p -> -1/(bare2+barp2)(
                    barp2p(Rv2-RT2)-4 Pi G a^2/H(bare1+barp1)(bare2+barp2)(Rv1-Rv2)+
                    + 2H/3 nabla2 Pi2 )
     }

  zero[2]=Collect[zero[1]//.subst1,{Y,Yp,Ypp,Rv1,Rv1p,RT1,RT1p,Rv2,Rv2p,RT2,RT2p},FullSimplify]

(* insert 9.51 *)

  testzero = (bare1+barp1)(
	     -3 RT1p - nabla2 Rv1/H + 12 Pi G a^2/H (
	     (bare1+barp1)(Rv1-RT1) + (bare2+barp2)(Rv2-RT2)  ))

  zero[3]=Collect[ (zero[2]-testzero)/.subst1,{Rv1,Rv1p,Rv2,Rv2p,RT1,RT1p,RT2,RT2p},FullSimplify]

  Print[zero[3]]

(**************************************************)


