function [ver,cst1,BC,BC_,icst,icst_] = verfun(SS,XX,CC,cstv,csts,cstn,bndv,bndn,i1,i2,par) %#ok<INUSL>

e  = SS(1);

z  = SS(2);

if i1==1
    bnd1 = '0';
elseif i1==par.n1
    bnd1 = '2';
else
    bnd1 = '1';
end

if i2==1
    bnd2 = '0';
elseif i2==par.n2
    bnd2 = '2';
else
    bnd2 = '1';
end

bnd = [bnd1,bnd2];

ver = base2dec(bnd,par.base)+1;

cst1  = NaN(par.ncc,1);

icst  = NaN(par.ncc,1);
icst_ = NaN(par.ncc,1);
BC    = NaN(par.ncc,1);
BC_   = NaN(par.ncc,1);


