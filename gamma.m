(* ::Package:: *)

(* ::Input:: *)
(*\[CapitalGamma]/:\[CapitalGamma][\[Omega]1_,A1_]\[CapitalGamma][\[Omega]2_,A2_]:=\[CapitalGamma][\[Omega]1 \[Omega]2,A1 +A2]*)
(*\[CapitalGamma]/:\[CapitalGamma][\[Omega]1_,A1_]\[Congruent]\[CapitalGamma][\[Omega]2_,A2_]:=Simplify[\[Omega]1 -\[Omega]2==0&&A1- A2==0]*)
(*\[CapitalGamma]Collect[\[CapitalGamma][\[Omega]_,A_]]:=\[CapitalGamma][Together[\[Omega]]//Expand,Collect[A,Subscript[r, _],Collect[#,Subscript[c, _],Together]&]]*)
(*\[CapitalGamma]Format[\[CapitalGamma][\[Omega]_,A_]]:=Module[{S,M},*)
(*S=Union@Cases[\[CapitalGamma][\[Omega],A],Subscript[(r|c), a_]:> a,Infinity];*)
(*M=Outer[Factor[\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \( *)
(*\*SubscriptBox[\(c\), \(#1\)], *)
(*\*SubscriptBox[\(r\), \(#2\)]\)]A\)]&,S,S];*)
(*M=Prepend[M,Subscript[r, #]&/@S]//Transpose;*)
(*M=Prepend[M,Prepend[Subscript[c, #]&/@S,\[Omega]]];*)
(*M//MatrixForm*)
(*]*)
(*ShiftT[P_]:=If[P===0,P,P*t^Exponent[P,1/t]//Expand] (*Shift the power of t to make the lowest power of t to be (t^0).*)*)
(*\[CapitalGamma]Format2[\[CapitalGamma][\[Omega]_,A_]]:=Module[{S,M},*)
(*S=Union@Cases[\[CapitalGamma][\[Omega],A],Subscript[(r|c), a_]:> a,Infinity];*)
(*M=Outer[ShiftT@Expand@Factor[ \!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \( *)
(*\*SubscriptBox[\(c\), \(#1\)], *)
(*\*SubscriptBox[\(r\), \(#2\)]\)]A\)]&,S,S];*)
(*M=Prepend[M,Subscript[r, #]&/@S]//Transpose;*)
(*M=Prepend[M,Prepend[Subscript[c, #]&/@S,ShiftT@\[Omega]]];*)
(*M//MatrixForm*)
(*]*)
(*Subscript[X, a_,b_]:=\[CapitalGamma][1,Subscript[r, a] Subscript[c, a]+(1-t)Subscript[r, a] Subscript[c, b]+t Subscript[r, b] Subscript[c, b]]*)
(*Subscript[\!\(\*OverscriptBox[\(X\), \(_\)]\), a_,b_]:=Subscript[X, a,b]/.t->t^-1*)
(*Subscript[m, a_,b_->k_][\[CapitalGamma][\[Omega]_,A_]]:=\[CapitalGamma][\[Omega](1-\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), *)
(*SubscriptBox[\(r\), \(a\)]]\( *)
(*\*SubscriptBox[\(\[PartialD]\), *)
(*SubscriptBox[\(c\), \(b\)]]A\)\)),A+((\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), *)
(*SubscriptBox[\(r\), \(a\)]]A\))(\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), *)
(*SubscriptBox[\(c\), \(b\)]]A\)))/(1-\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), *)
(*SubscriptBox[\(r\), \(a\)]]\( *)
(*\*SubscriptBox[\(\[PartialD]\), *)
(*SubscriptBox[\(c\), \(b\)]]A\)\))]/.{(Subscript[c, b]|Subscript[r, a])->0,Subscript[c, a]->Subscript[c, k],Subscript[r, b]->Subscript[r, k]}//\[CapitalGamma]Collect(*//\[CapitalGamma]Collect*)*)
(*Subscript[II, i_]:=\[CapitalGamma][1,Subscript[r, i] Subscript[c, i]+Subscript[r, -i] Subscript[c, -i]]*)
(**)
(*Subscript[XXp, i_,j_]:=\[CapitalGamma][1,(1-t) Subscript[c, -j] Subscript[r, -i]+(Subscript[c, -i]+(1-t) Subscript[c, j]) Subscript[r, -i]+Subscript[c, i] Subscript[r, i]+t Subscript[c, -j] ((1-t) Subscript[r, i]+t Subscript[r, -j])+t Subscript[c, j] ((1-t) Subscript[r, i]+t Subscript[r, j])]*)
(*Subscript[XXn, i_,j_]:=\[CapitalGamma][1,Subscript[c, -i] Subscript[r, -i]+(Subscript[c, i]+(1-1/t) Subscript[c, -j]) Subscript[r, i]+(1-1/t) Subscript[c, j] Subscript[r, i]+(Subscript[c, -j] ((1-1/t) Subscript[r, -i]+Subscript[r, -j]/t))/t+(Subscript[c, j] ((1-1/t) Subscript[r, -i]+Subscript[r, j]/t))/t]*)
(*Subscript[mm, i_, j_->k_][G_]:=G//Subscript[m, i,j->k]//Subscript[m, -i,-j->-k]//\[CapitalGamma]Collect*)
(*Subscript[HHhv, a_,b_][h_]:=\[CapitalGamma][1,t Subscript[c, -a] Subscript[r, -a]+(1-t) Subscript[c, -b] Subscript[r, -a]+(t-h t) Subscript[c, -b] Subscript[r, a]+h Subscript[c, b] Subscript[r, a]+(1-t) Subscript[c, -a] Subscript[r, -b]+Subscript[c, a] Subscript[r, -b]+(-t+h t) Subscript[c, -b] Subscript[r, -b]+(1-h) Subscript[c, b] Subscript[r, -b]+t Subscript[c, -b] Subscript[r, b]]*)
(**)
(*(*"merge with" operations for repeatedly merging in a sequence into one large strand with label s*)*)
(*Subscript[mXp, i_,j_,k_][G_]:= G Subscript[XXp, i,j]// Subscript[mm, s, k->s]*)
(*Subscript[mXn, i_,j_,k_][G_]:= G Subscript[XXn, i,j]// Subscript[mm, s, k->s]*)
(*Subscript[mH, h_,i_,j_,k_][G_]:= G Subscript[HHhv, i,j][h]//Subscript[mm, s, k->s]*)
(**)
