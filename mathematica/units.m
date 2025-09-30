
(*unit conversions*)
 
  mpl = 1.2209 10^(19) GeV;
  c = 299792458 m/s; 
  hbar = 6.582119569 10^(-25) GeV s;
  pc = 3.0857 10^(16) m;
  ly = 0.946073 10^(16) m;
  K = 8.617333 10^(-5) 10^(-9) GeV;
  Hz = 1/s;
  H0oh = 10^5 m/s /(10^6 pc);

(*numerical values*)

  T0 = 2.7255 K;
  s0T3 = 3.92 2 Pi^2/45;
  Tl = 10^(-3) GeV;
  elT4 = 3.4945;
  slT3 = 4.6565;
  Te = 10^(15) GeV;
  eeT4 = 34.602;
  seT3 = 46.135;
  Tswitch = 10^(-12) mpl;
  erT4 = 34.421;
  srT3 = 45.898;

(*hubble rate*)

  Print["h/H0/pc:"]
  Print[N[c/H0oh/pc]]
  Print["h/H0/ly:"]
  Print[N[c/H0oh/ly]]
  Print["H0/h/Hz:"]
  Print[N[H0oh/Hz]]
  Print[" "]

(*e-folds*)

  Print["mpl*Mpc:"]
  Print[N[mpl 10^6 pc /(c hbar)]]
  Print[" "]

  Print["e-folds from Te to T0:"]
  Print[N[Log[Te/T0 (seT3/s0T3)^(1/3)]]]
  Print[" "]

  Print["e-folds from Tswitch to T0:"]
  Print[N[Log[Tswitch/T0 (srT3/s0T3)^(1/3)]]]
  Print[" "]

(*redshift*)

  Print["redshifts:"]
  redshifte=(seT3/s0T3)^(1/3) Sqrt[3/(8Pi)] (mpl/Te)/Sqrt[eeT4];
  Print[N[redshifte]]
  redshiftl=(slT3/s0T3)^(1/3) Sqrt[3/(8Pi)] (mpl/Tl)/Sqrt[elT4];
  Print[N[redshiftl]]
  Print[" "]

(*frequency*)

  Print["test frequencies:"]
  Print[N[2 Pi hbar/s/T0]]
  Print[N[s T0/(2 Pi hbar)/redshifte]]
  Print[N[s T0/(2 Pi hbar)/redshiftl]]
  Print[" "]

(*wavelength*)

  Print["test wavelengths:"]
  Print[N[2 Pi hbar c/pc/T0]]
  Print[N[2 Pi hbar c/pc/T0 redshifte]]
  Print[N[2 Pi hbar c/pc/T0 redshiftl]]
  Print[" "]

(*relation between frequency and wavelength*)

  Print["relation between f_0 and lambda_0:"]
  Print[N[c s/pc]]
  Print[" "]

(*matter-radiation equality*)

  Print["matter-radiation equality:"]
  Print[N[(1+5.4) GeV 2 Zeta[3]/Pi^2 6.1 10^(-10)/(Pi^2/30 (2 +7/4(4/11)^(4/3) 3.044))]]
  Print[" "]
