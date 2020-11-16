if vtmp==1
    myfun = @(X) Ffun00(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun00(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun00(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==2
    myfun = @(X) Ffun01(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun01(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun01(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==3
    myfun = @(X) Ffun02(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun02(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun02(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==4
    myfun = @(X) Ffun10(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun10(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun10(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==5
    myfun = @(X) Ffun11(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun11(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun11(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==6
    myfun = @(X) Ffun12(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun12(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun12(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==7
    myfun = @(X) Ffun20(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun20(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun20(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==8
    myfun = @(X) Ffun21(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun21(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun21(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==9
    myfun = @(X) Ffun22(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun22(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun22(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==0
    myfun = @(X) Ffun33(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    dFF   = dFfun33(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
    FF    =  Ffun33(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
end

