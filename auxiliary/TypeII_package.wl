(* ::Package:: *)

(*Add path to the perturbative data which will be used to initialise*)

PerturbativeDataPath = "../data/perturbative/"

(*uncomment when TypeII_script.wls is run in manual mode*)

(*PerturbativeDataPath = NotebookDirectory[]<>"../data/perturbative/"*)


(*Oscillator numbers*)

nOsc = {{nb[1],nb[2]},{nf[1],nf[2],nf[3],nf[4]},{na[1],na[2]}}


(*Length*)
L := 1/2 (Sum[nf[j], {j, 1, 4}] + Sum[na[j], {j, 1, 2}] - Sum[nb[j], {j, 1, 2}]);

(*Bare dimension*)
\[CapitalDelta]0 := 1/2 Sum[nf[j], {j, 1, 4}] + Sum[na[j], {j, 1, 2}];


(*Compute asymptotics of Pa and Qi*)

\[Lambda]0[a_] := nf[a]+{2,1,0,-1}[[a]];
\[Lambda][a_] := \[Lambda]0[a] + \[CapitalLambda];
\[Nu]0[i_] := {-L - nb[1] - 1, -L - nb[2] - 2, na[1] + 1, na[2]}[[i]];
\[Nu][i_] := \[Nu]0[i] + (\[CapitalDelta] - \[CapitalDelta]0)/2 {-1,-1,1,1}[[i]]- \[CapitalLambda];

\[CapitalLambda] = (1 - \[Lambda]0[1] - \[Lambda]0[4])/2

powP = - (Table[\[Lambda][a], {a, 1, 4}])
powQ = - (Table[\[Nu][i], {i, 1, 4}]) - 1


(*Function that generates the main substitution rule sb0, which contains the current values of all the parameters \[CapitalDelta] and {c}*)
(*Format of the vector of parameters = {g,\[CapitalDelta],Subscript[c, 1,1],Subscript[c, 1,2],...,Subscript[c, 1,cutP],Subscript[c, 2,1],Subscript[c, 2,2],...,Subscript[c, 2,cutP],Subscript[c, 3,1],Subscript[c, 3,2],...,Subscript[c, 3,cutP],Subscript[c, 4,1],Subscript[c, 4,2],...,Subscript[c, 4,cutP],}*)

\[DoubleStruckCapitalF]V0[Plist_List]:=Block[{},
If[Length[Flatten[Plist]]!=4cutP+2,Print["Length error"];Return[];];
sb0=({
c[2,-(powP[[1]]-powP[[2]])]->0,
c[3,-(powP[[1]]-powP[[3]])]->0,
c[3,-(powP[[2]]-powP[[3]])]->0,
c[4,-(powP[[1]]-powP[[4]])]->0,
g->Plist[[1]]//N[#,WP]&,
\[CapitalDelta]->Plist[[2]]//N[#,WP]&,
Sequence@@Table[c[1,i]->Plist[[2 + i]]//N[#,WP]&,{i,1,cutP,1}],
Sequence@@Table[c[2,i]->Plist[[2 + cutP + i]]//N[#,WP]&,{i,1,cutP,1}],
Sequence@@Table[c[3,i]->Plist[[2 + cutP +cutP + i]]//N[#,WP]&,{i,1,cutP,1}],
Sequence@@Table[c[4,i]->Plist[[2 + cutP + cutP + cutP + i]]//N[#,WP]&,{i,1,cutP,1}]
}//SetPrecision[#,WP]&)/.c[aa_,bb_]:>c[Rationalize[aa],Rationalize[bb]]//Association;
]


(*Find fixed and independent paramaters by removing those parameters which are set to zero by gauge fixing and the parameter fixed by the constraint*)

IndepParams={\[CapitalDelta],Table[c[1,i],{i,1,cutP,1}],Table[c[2,i],{i,1,cutP,1}],Table[c[3,i],{i,1,cutP,1}],Table[c[4,i],{i,1,cutP,1}]}/.c[2,-(powP[[1]]-powP[[2]])]->0/.c[3,-(powP[[1]]-powP[[3]])]->0/.c[3,-(powP[[2]]-powP[[3]])]->0/.c[4,-(powP[[1]]-powP[[4]])]->0/.c[1,na[1]-na[2]+1]->0//Flatten//DeleteCases[#,0]&;
ParamsInPlist={g,\[CapitalDelta],Table[c[1,i],{i,1,cutP,1}],Table[c[2,i],{i,1,cutP,1}],Table[c[3,i],{i,1,cutP,1}],Table[c[4,i],{i,1,cutP,1}]}//Flatten;
FixedParams=Complement[ParamsInPlist,IndepParams];
nnns=Table[Position[ParamsInPlist,ip],{ip,IndepParams}]//Flatten;
nsFixed=Table[Position[ParamsInPlist,ip],{ip,FixedParams}]//Flatten//Sort;


(*Adjust the size of the vector of parameters to be 4 cutP + 2*)

IncrBs[Plist_List]:=Block[{},
If[Length[Flatten[Plist]]==4cutP+2,
coefs = Plist;
Return[coefs]
];
If[Length[Flatten[Plist]]!=4cutP+2,
lst=Drop[Plist,2];
oldcutP = Length[lst]/4;
lstc12 = Take[lst,  2 oldcutP];
lstc34 = Take[lst,- 2 oldcutP];
lstc1=Take[lstc12,oldcutP];
lstc2 =Take[lstc12,-oldcutP];
lstc3=Take[lstc34,oldcutP];
lstc4=Take[lstc34,-oldcutP];
lstc1=Take[lstc1~Join~Table[0,{cutP - oldcutP}],cutP];
lstc2=Take[lstc2~Join~Table[0,{cutP - oldcutP}],cutP];
lstc3=Take[lstc3~Join~Table[0,{cutP - oldcutP}],cutP];
lstc4=Take[lstc4~Join~Table[0,{cutP - oldcutP}],cutP];
Return[Take[Plist,2]~Join~lstc1~Join~lstc2~Join~lstc3~Join~lstc4]
]
]


(*Function to boost precision*)

GoodB[Plist_List]:=Plist// SetPrecision[#,WP*3/2]&;


(*Function to interpolate parameters using Fit to produce a vector of starting points at g = g0*)
(*The function uses eOrder previous saved points which have at least 10^(-prec) precision*)

InterpolateIn[g0_,eOrder_,prec_]:=Block[{},
val=g0;
saved0=Select[saved,error[#]<10^(-prec)&]//Union;
extOrder=Min[{eOrder,Length[saved0]}];
savedNearest=Take[SortBy[saved0,Abs[#-val]^2&],extOrder];
savedNearest=savedNearest//DeleteDuplicates;
Bss=Table[IncrBs[params[sn]],{sn,savedNearest}];
fit=(Fit[Table[{N[savedNearest[[a]],2WP],bb[a]},{a,extOrder}],x^Range[0,Max[0,extOrder-1]],x]/.x->val);
Bs0=GoodB[fit/.bb[a_]:>Bss[[a]]];
Bs0[[1]]=val;
Bs0
];


(*Generate the initial vector of parameters from the perturbative data of arXiv:1812.09238*)
(*In practice, we use the substitution rule sbWeak which is saved in the perturbative data files on GitHub*)
(*Below we assume that the perturbative data file is located in a directory called data/perturbative and the directory data/perturbative is located in the parent directory of this file*)
(*The path at fnameA would need to be adjusted accordingly should the reader choose a different location for the perturbative data*)

FromPert[gloop_]:=Block[{res},
fnameA = PerturbativeDataPath<>"perturbative_data_Delta0"<>ToString[\[CapitalDelta]0]<>"_b1"<>ToString[nb[1]]<>"_b2"<>ToString[nb[2]]<>"_f1"<>ToString[nf[1]]<>"_f2"<>ToString[nf[2]]<>"_f3"<>ToString[nf[3]]<>"_f4"<>ToString[nf[4]]<>"_a1"<>ToString[na[1]]<>"_a2"<>ToString[na[2]]<>"_sol"<>ToString[sol]<>".mx";
	nexpn = 50;
	Clear[cg,\[CapitalDelta]g];
	\[CapitalDelta]Anz = Sum[\[CapitalDelta]g[i]g^i,{i,0,nexpn,2}];
	\[CapitalDelta]g[0]=\[CapitalDelta]0;
    cg[2,-(powP[[1]]-powP[[2]])]=0;
    cg[3,-(powP[[1]]-powP[[3]])]=0;
    cg[3,-(powP[[2]]-powP[[3]])]=0;
    cg[4,-(powP[[1]]-powP[[4]])]=0;
	cg[1,k_]:= Sum[cg[1,k,n]g^n,{n,-6,nexpn,1}]//Collect[#,g]&;
	cg[2,k_]:= Sum[cg[2,k,n]g^n,{n,-6,nexpn,1}]//Collect[#,g]&;
	cg[3,k_]:= Sum[cg[3,k,n]g^n,{n,-6,nexpn,1}]//Collect[#,g]&;
	cg[4,k_]:= Sum[cg[4,k,n ]g^n,{n,-6,nexpn,1}]//Collect[#,g]&;
	Get[fnameA];
	\[DoubleStruckCapitalF]V0[invsb00=a/@Range[4 cutP+2]];
	invsb0=invsb00/. (Rationalize[Normal[sb0]]/. (a_->b_):>b->a);
	res=GoodB[IncrBs[invsb0/. \[CapitalDelta]->\[CapitalDelta]Anz/. c->cg/.sbWeak/.\[CapitalDelta]g[_]->0/.cg[___]->0/.g->gloop]];
res[[1]]=gloop;
res
]


(*Function to generate starting points at weak coupling using oa combination of perturbative data, and previous saved numerical data*)
(*The function uses eOrder previous saved points which have at least 10^(-prec) precision*)

InterpolateWeak[g0_,eOrder_,prec_]:=Block[{ord=eOrder},
val=g0;
saved0=Union[Select[saved,error[#1]<10^-prec&]];
extOrder=Min[{eOrder,Length[saved0]}];
savedNearest=Take[SortBy[saved0,Abs[#1-val]^2&],extOrder];
savedNearest=DeleteDuplicates[savedNearest];
Bss=Table[IncrBs[params[sn]],{sn,savedNearest}];
fw=FromPert[g];
\[DoubleStruckCapitalF]V0[invsb00=a/@Range[4 cutP+2]];
invsb0=invsb00/. (Rationalize[Normal[sb0]]/. (a_->b_):>b->a);
forfit=invsb0/. \[CapitalDelta]->\[CapitalDelta]Anz/. c->cg/.sbWeak/. \[CapitalDelta]g->cg;
If[extOrder==0,
Return[fw/. g->gloop]
];
Bs0=Join[
{g0},
Table[
If[!MemberQ[nnns,trm],
0,
vars=Take[Union[Cases[forfit[[trm]],cg[__],\[Infinity]]],extOrder];
tovar=Apply[Rule,Transpose[{vars,Table[va[i],{i,extOrder}]}],{1}];
ans=forfit[[trm]]/. tovar/. cg[__]->0;
findfit=FindFit[Transpose[{Transpose[Bss][[1]],Transpose[Bss][[trm]]}],ans,Union[Join[Cases[{ans},\[CapitalDelta]g[_],\[Infinity]],Cases[{ans},va[__],\[Infinity]]]],g];
Expand[ans/. findfit/. g->g0]]
,{trm,2,Length[fw]}]
];
Bs0[[1]]=val;
res=GoodB[IncrBs[Bs0]];
res[[1]]=gloop;
res
];
