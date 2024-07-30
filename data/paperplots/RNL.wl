(* ::Package:: *)

(* RNL *)

BeginPackage["RNL`"]

fst::usage = "fst[a] = a[[1]]"
snd::usage = "snd[a] = a[[2]]"
Zip::usage = "Zip[a,b,...] = Transpose[{a,b,...}]"
Eqlen::usage = "Truncates all lists in list to length of shortest list"
ilist::usage = "ilist[l] = Range[Length[l]]"
fn::usage = "fn[r] = Replace[#,r] &"
NumberFromString::usage = "Outputs number corresponding to string"
LinSpace::usage = "Numbers equally spaced linearly"
LogSpace::usage = "Numbers equally spaced logarithmically"
trapezoid::usage = "trapezoid[{{x_i,y_i}}] = trapezoid rule integration of points"
FindRootBisect::usage = "FineRootBisect[f,x0,x1,n] finds root of f[x] in interval [x0,x1] by n rounds of bisection"
applyN::usage = "applyN[f,a] = {f[1,a], f[2,a], ...}"
FindMinGolden::usage = "FindMinGolden[f,x0,x1,n] finds minimum of f[x] in bracketing interval [x0,x1] by n rounds of Golden Ratio search"
importCSV::usage = "ImportString[s,\"CSV\"]"

(* ************************************************************** *)

Begin["`Private`"]

fst[a_] := a[[1]]; snd[a_] := a[[2]];
Zip[a___] := Transpose[{a}]
Eqlen[l_] := With[{lmin = Min[Length /@ l]}, Take[#, lmin] & /@ l]
ilist[l_] := Range[Length[l]]
fn[rule_] := Replace[#, rule] &

(* TODO - there's almost certainly a better way of doing this,
but going with this ftm ... *)
NumberFromString[s_] := Read[StringToStream[s],Number]

(* This probably doesn't quite belong here, but put it here ftm *)

LinSpace[a_,b_,n_] := Range[a,b,(b-a)/(n-1)]
LogSpace[a_,b_,n_] := With[{la = Log[a], lb = Log[b]},Exp[Range[la,lb,(lb-la)/(n-1)]]]

(* Should do a better quadrature rule here ... *)
trapezoid[t_] := 1/2 Total[Differences[t[[All,1]]] ListCorrelate[{1,1},t[[All,2]]]]

FindRootBisect[f_, x0_, x1_, n_] := Module[{xa1, xb1, ftmp},
  {xa1, xb1} = If [f[x0] < 0, {x0, x1}, {x1, x0}];
  Assert[f[xa1] <= 0 && f[xb1] >= 0];
  ftmp[xa_, xb_, m_] := With[{xc = 0.5 (xa + xb)}, If[m == 1, {xa, xb},
     If[f[xc] <= 0, ftmp[xc, xb, m - 1], ftmp[xa, xc, m - 1]]]];
  ftmp[xa1, xb1, n]
  ]

FindMinGolden[f_,x0_,x1_,n_] := Module[{phi,ftmp,x2},
phi = N[GoldenRatio];
x2 = x0 + (x1-x0)/phi;
ftmp[xa_,xb_,xc_,m_] := If[m == 1, {xa,xc},
If[(xb - xa) > (xc - xb),
With[{xn = xb - (xb-xa)/phi},
If[f[xn] < f[xb],
ftmp[xa,xn,xb,m-1], ftmp[xn,xb,xc,m-1]]],
With[{xn = xb + (xc-xb)/phi},
If[f[xn] < f[xb],
ftmp[xb,xn,xc,m-1],ftmp[xa,xb,xn,m-1]]]]];
ftmp[x0,x2,x1,n]]

applyN[f_, a_] := f /@ Zip[ilist[a], a];

importCSV[s_] := ImportString[s,"CSV"];

(* *)

End[]

(* ************************************************************** *)

EndPackage[]















