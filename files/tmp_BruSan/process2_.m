if vtmp==1
    L(:,i1,i2) = Lfun00(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==2
    L(:,i1,i2) = Lfun01(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==3
    L(:,i1,i2) = Lfun02(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==4
    L(:,i1,i2) = Lfun10(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==5
    L(:,i1,i2) = Lfun11(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==6
    L(:,i1,i2) = Lfun12(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==7
    L(:,i1,i2) = Lfun20(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==8
    L(:,i1,i2) = Lfun21(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==9
    L(:,i1,i2) = Lfun22(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);

elseif vtmp==0
    L(:,i1,i2) = Lfun33(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);
end

