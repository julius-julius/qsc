#!/usr/bin/env wolframscript
(* ::Package:: *)

(*Helper functions*)

AT = AbsoluteTime[];
SetOptions[$Output, FormatType->OutputForm];
ToString1=ToString[#,PageWidth->Infinity]&;


(*Path to TypeII_exec.out*)
ExecPath = "../core/";

(*Path to TypeII_package.wl*)
PackagePath = "";

(*Path to where you want to store the output*)
OutputPath = "../data/output/";

(*Import hyperparameters from TypeII_run.ipynb*)

Print["parameters in WL:",$ScriptCommandLine];

nb[1] = ToExpression[$ScriptCommandLine[[2]]];
nb[2] = ToExpression[$ScriptCommandLine[[3]]];
nf[1] = ToExpression[$ScriptCommandLine[[4]]];
nf[2] = ToExpression[$ScriptCommandLine[[5]]];
nf[3] = ToExpression[$ScriptCommandLine[[6]]];
nf[4] = ToExpression[$ScriptCommandLine[[7]]];
na[1] = ToExpression[$ScriptCommandLine[[8]]];
na[2] = ToExpression[$ScriptCommandLine[[9]]];
sol = ToExpression[$ScriptCommandLine[[10]]];

cutP = ToExpression[$ScriptCommandLine[[11]]];
nPoints = ToExpression[$ScriptCommandLine[[12]]];
cutQai = ToExpression[$ScriptCommandLine[[13]]];
QaiShift = ToExpression[$ScriptCommandLine[[14]]];
\[Epsilon]=10^(-ToExpression[$ScriptCommandLine[[15]]]);
WP=ToExpression[$ScriptCommandLine[[16]]];

gloop = ToExpression[$ScriptCommandLine[[17]]];

interOrder = ToExpression[$ScriptCommandLine[[18]]];
interMethod = ToExpression[$ScriptCommandLine[[19]]];
extrations=ToExpression[$ScriptCommandLine[[20]]];
If[extrations==1,extrations=5];
errorTol = ToExpression[$ScriptCommandLine[[21]]];

(*Import TypeII_package.wl*)
Get[PackagePath<>"TypeII_package.wl"]


(*(*uncomment to run in manual mode - to run without Python launcher*)


(*Path to TypeII_exec.out*)
SetDirectory[NotebookDirectory[]<>"../core/"];
ExecPath = "";

(*Path to TypeII_package.wl*)
PackagePath = NotebookDirectory[];

(*Path to where you want to store the output*)
OutputPath = NotebookDirectory[]<>"../data/output/";


nb[1]=0;
nb[2]=1;
nf[1]=2;
nf[2]=2;
nf[3]=1;
nf[4]=1;
na[1]=1;
na[2]=0;
sol = 1;

cutP = 16;
nPoints = cutP + 2;
cutQai = 30;
QaiShift = 100;
\[Epsilon] = 10^(-30);
WP = 186;

gloop = 1/1000;

interOrder = 10;
interMethod = 1;
kernels = 1;
errorTol = 36;

extrations=5;

Get[PackagePath<>"TypeII_package.wl"]*)


(*Import numerical data file if it exists*)
Print["Loading data..."];
fnameN= OutputPath<>"numerical_data_Delta0"<>ToString1[\[CapitalDelta]0]<>"_b1"<>ToString1[nb[1]]<>"_b2"<>ToString1[nb[2]]<>"_f1"<>ToString1[nf[1]]<>"_f2"<>ToString1[nf[2]]<>"_f3"<>ToString1[nf[3]]<>"_f4"<>ToString1[nf[4]]<>"_a1"<>ToString1[na[1]]<>"_a2"<>ToString1[na[2]]<>"_sol"<>ToString1[sol]<>".mx";

If[!FileExistsQ[fnameN],
saved={};
Clear[init,params,error,values,digits];,
Get[fnameN];
Print["Loaded:",fnameN];
];


(*generate starting points via extrapolation*)

Print["Extrapolation..."];
If[
	interMethod===1,
	Print["using::perturbative"];
	BS[0]= InterpolateWeak[gloop,interOrder,20];,
	Print["using::extrapolation"];
	BS[0]=InterpolateIn[gloop,interOrder,20];
]
Print["gmax=",saved//N//Re//Max];


(*conversion of hyperparameters to internal format of C++ code*)

init11=SetPrecision[BS[0],WP];

g0=init11[[1]];
\[Gamma]0=init11[[2]]-\[CapitalDelta]0;(* only the anomalous part of \[CapitalDelta]*)
Mt = Table[\[Lambda][a],{a,4}];

c1=1/g0^Mt[[1]]*Table[init11[[2+j]],{j,1,cutP}];
c2=1/g0^Mt[[2]]*Table[init11[[2+cutP+j]],{j,1,cutP}];
c3=1/g0^Mt[[3]]*Table[init11[[2+2*cutP+j]],{j,1,cutP}];
c4=1/g0^Mt[[4]]*Table[init11[[2+3*cutP+j]],{j,1,cutP}];
c={Insert[c1,1,1],Insert[c2,1,1],Insert[c3,1,1],Insert[c4,1,1]};
(* A false  value (=1) for c_{a,0} is inserted, which is irrelevant for the C++ program, since there c_{a,0} is computed from the input value of \[CapitalDelta] *)

convertionFactor = Table[1/g0^Mt[[a]],{a,4},cutP]//Flatten;

precssf=2*(errorTol-4); (* iteration qiuting condition parameter,*)
maxiter=extrations;  (* iteration qiuting condition parameter, if the number of iterations is larger than maxiter, than the programm finishes running, and return the value, it obtained upto maxiter iterations*)


(*generate the string to attach to the executable file TypeI_exec.out*)

ingen=Table["",{k,1,24+4*cutP}];
For[k=0,k<=cutP,k++,     
ingen[[21+k]]=ToString1[CForm[Im[c[[1,k+1]]]]];
ingen[[22+cutP+k]]=ToString1[CForm[Re[c[[2,k+1]]]]];
ingen[[23+2*cutP+k]]=ToString1[CForm[Im[c[[3,k+1]]]]];
ingen[[24+3*cutP+k]]=ToString1[CForm[Re[c[[4,k+1]]]]];
];  
(**********************************************)
(**********************************************)
ingen[[1]]=ToString1[WP];(*working precision*)
ingen[[2]]=ToString1[cutP];
ingen[[3]]=ToString1[cutQai];
ingen[[4]]=ToString1[QaiShift];
ingen[[5]]=ToString1[nPoints];
ingen[[6]]=ToString1[-Round[Log10[\[Epsilon]]]];(*exponent of \[Epsilon]*)
ingen[[7]]=ToString1[precssf];(*number of iterations*)
ingen[[8]]=ToString1[errorTol];(*last change in Delta*)
ingen[[9]]=ToString1[maxiter];(*max iterations*)
ingen[[10]]=ToString1[nb[1]];
ingen[[11]]=ToString1[nb[2]];
ingen[[12]]=ToString1[nf[1]];
ingen[[13]]=ToString1[nf[2]];
ingen[[14]]=ToString1[nf[3]];
ingen[[15]]=ToString1[nf[4]];
ingen[[16]]=ToString1[na[1]];
ingen[[17]]=ToString1[na[2]];
ingen[[18]]=ToString1[sol];(*multiplicity label sol*)
ingen[[19]]=ToString1[NumberForm[Re[N[g0,{WP,WP}]],{WP,WP}]];
ingen[[20]]=ToString1[CForm[Re[N[\[Gamma]0,{WP,WP}]]]];


(*RUNNING C++*)


(*auxiliary function*)
Fprec[x_,prec_]:=N[Floor[x*10^(15-(MantissaExponent[x][[2]]-1))]/10^15,{prec,prec}]*10^(MantissaExponent[x][[2]]-1);


Err:=(Print["ERROR: no output created"];
Print["PRECISION||",11//Log10//N,"||"];
Print["EXTRA||",11//Log10//N,"||"];
Print["icurrent||",440,"||"];
Print["ATIME||", Round[10*(AbsoluteTime[]-AT)]*0.1,"||"];
Quit[];);


torun="! ./"<>ExecPath<>"TypeII_exec.out";
(**)
For[k=1,k<=Length[ingen],k++,torun=torun<>" "<>ingen[[k]] ];


(*Print[torun]*)


Runingen[ingen_]:=Block[{},
Print["Running C++..."];
output=Import[torun,"String"];
Toexpr=ToExpression[output];
Toexpr[[2]]=Fprec[Toexpr[[2]],WP];
Toexpr
]


Runingen[ingen];


If[Length[Toexpr]<8,Err];


(*OUTPUT CONVERSION*)


(* ::Input:: *)
(*(* the output structure is as follows *)*)


(* ::Input:: *)
(*(* Toexpr[[1]] is an errorcode: if it =0 then the result is within the prescribed precisions , if it is =1, maxiter is reached, and the required precisons are not reached *)*)


(* ::Input:: *)
(*(* anomalous  part of \[CapitalDelta] *)*)


\[Gamma]=Toexpr[[4]];


(* ::Input:: *)
(*(* c_{1,n} n=0...,cutP/2*)*)


c1=Drop[Toexpr[[5]],1];
c2=Drop[Toexpr[[6]],1];
c3=Drop[Toexpr[[7]],1];
c4=Drop[Toexpr[[8]],1];


BSicurrent={{gloop,\[CapitalDelta]0+\[Gamma]},Flatten[{c1,c2,c3,c4}]/convertionFactor}//Flatten;


(* ::Input:: *)
(*(* further output and running parameter values *)*)


others=Toexpr[[9]];
vec1=Toexpr[[10]];
vec2=Toexpr[[11]];
Print["Stest||",N[Log10[#//Abs]]&/@vec2,"||"];
Print["Shifttext||",N[Log10[#//Abs]]&/@vec1,"||"];


(* ::Input:: *)
(*(* number of iterations done *)*)


icurrent=others[[1]];


(* ::Input:: *)
(*(*   "the sum of |\[Delta]c_{a,n}|^2 *)*)


(* ::Input:: *)
(*(*  |\[CapitalDelta]^(n+1)-\[CapitalDelta]^(n)| at the end of running *) *)


d\[CapitalDelta]=others[[2]];


(*not ideal - needs improvement*)
\[DoubleStruckCapitalS]best=If[Toexpr[[1]]===0,Abs[d\[CapitalDelta]],2];
Print["Log10[Sbest]=",\[DoubleStruckCapitalS]best//N//Log10];


(*OUTPUT*)


combinedPrec = 10^Max[{\[DoubleStruckCapitalS]best//Log10,Abs[others[[13]]]//Log10}];
combinedPrec = Last[(#//Abs)&/@vec2]^2;
\[DoubleStruckCapitalS]best=combinedPrec;


BS[icurrent]=BSicurrent;


DelOld=BS[0][[2]];
DelNew=BS[icurrent][[2]];
prec=Abs[DelNew-DelOld+10^(-WP+10)];


valdata={{cutP,nPoints},{cutQai,QaiShift},{Log[10,\[Epsilon]],WP}};

If[!FileExistsQ[fnameN],
saved={};
Clear[params,error,values,digits];
];

PW0=WP;

If[combinedPrec<10^(-errorTol)&&prec<1/100,
Quiet[Get[fnameN];];
toadd=gloop;
Print["Saved: ",toadd//N,",  \[CapitalDelta]=",BS[icurrent][[2]]//N//Re,",  \[DoubleStruckCapitalS]=",(error[toadd]=\[DoubleStruckCapitalS]best//N)//Log10];
init[toadd]=N[BS[0],PW0];
params[toadd]=N[BS[icurrent],PW0];
values[toadd]=valdata;
digits[toadd]=prec;
command[toadd]=$ScriptCommandLine;
saved=Union[Append[saved,toadd]];
Quiet[DeleteFile[fnameN];];
Save[fnameN,{saved,init,params,error,values,digits,command}];
savedQ = 1;
,
savedQ = 0]


Print["PRECISION||",\[DoubleStruckCapitalS]best//Log10//N,"||"]
Print["EXTRA||",Abs[others[[13]]]//Log10//N,"||"]
Print["iniS||",First[N[Log10[#//Abs]]&/@vec2],"||"]
Print["resS||",Last[N[Log10[#//Abs]]&/@vec2],"||"]
Print["dDelta||",N[Log10@prec],"||"]
Print["resDelta||",N[DelNew//Re],"||"]
Print["itr||",icurrent,"||"]
Print["shifterror||",Last[N[Log10[#//Abs]]&/@vec1],"||"]
Print["newtonGoal||",N[Log10[others[[2]]]],"||"]
Print["flagError||",Toexpr[[1]],"||"]
Print["flagSaved||",savedQ,"||"]
Print["atime||", Round[100000*(AbsoluteTime[]-AT)]*0.00001,"||"];


Print["end wls"]
