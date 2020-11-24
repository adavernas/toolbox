function [dFvars,Fvars,Xvars,Lvars,CCvars] = writefun(param,state,delta,value,...
                    vars,vars_,last,varl0,varm1,varm2,varp1,varp2,bnd,cstv,par)
 
%% Parameters
[~,npr] = size(param);
for ip=1:npr
    eval([param{ip},' = par.',param{ip},';'])
end
 
ns  = par.ns;
nd  = par.nd;
nv  = par.nv;
nx  = par.nx;
nx_ = par.nx_;
nl  = par.nl;
ncc = par.ncc;
np  = par.np;
 
if par.ns==1
all = [param,state,delta,value,vars,vars_,last,varl0,varm1,varp1];
for i=1:np+ns+nd+nv+nx+nx_+4*nl
    eval([all{i},' = sym(''',all{i},''');'])
end
elseif par.ns==2
all = [param,state,delta,value,vars,vars_,last,varl0,varm1,varm2,varp1,varp2];
for i=1:np+ns+nd+nv+nx+nx_+6*nl
    eval([all{i},' = sym(''',all{i},''');'])
end
end
 
Vzero = sym('Vzero');
Vone  = sym('Vone');
 
%% System of Equations for q
if str2double(bnd(1))==1 
    
    qe   = ( qpe*dme   - q*dme   + q*dpe   - qme*dpe   )/(2*dpe*dme);
    qee  = ( qepe*dme  - qe*dme  + qe*dpe  - qeme*dpe  )/(2*dpe*dme);
    qeee = ( qeepe*dme - qee*dme + qee*dpe - qeeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    qe   = qepe - qeepe*dpe;
    qee  = qeepe - qeeepe*dpe;
    qeee = ( qeepe - qee )/dpe;
 
elseif str2double(bnd(1))==2
    qe   = qeme + qeeme*dme;
    qee  = qeeme + qeeeme*dme;
    qeee = ( qee - qeeme )/dme;

end
 
if str2double(bnd(2))==1       
    qz   = ( qpz*dmz   - q*dmz   + q*dpz   - qmz*dpz  )/(2*dpz*dmz);
    qzz  = ( qzpz*dmz  - qz*dmz  + qz*dpz  - qzmz*dpz )/(2*dpz*dmz);
    qzzz = ( qzzpz*dmz - qzz*dmz + qzz*dpz - qzzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    qz   = qzpz - qzzpz*dpz;
    qzz  = qzzpz - qzzzpz*dpz;
    qzzz = ( qzzpz - qzz )/dpz;
    
elseif str2double(bnd(2))==2
    qz   = qzmz  + qzzmz*dmz;
    qzz  = qzzmz + qzzzmz*dmz;
    qzzz = ( qzz - qzzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    qze_ = ( qzpe*dme - qz*dme + qz*dpe - qzme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    qze_ = qezpe - qezepe*dpe;
    
elseif str2double(bnd(1))==2
    qze_ = qezme + qezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    qez_ = ( qepz*dmz - qe*dmz + qe*dpz - qemz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    qez_  = qezpz - qezzpz*dpz;
 
elseif str2double(bnd(2))==2
    qez_  = qezmz + qezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    qe   = 0;
    qee  = 0;
    qeee = 0;
    
    qz   = 0;
    qzz  = 0;
    qzzz = 0;
    
    qez  = 0;
    qeze = 0;
    qezz = 0;
 
else
    qez = 1/2*qez_ + 1/2*qze_;
end
 
if str2double(bnd(1))==1
    qeze  = ( qezpe*dme  - qez*dme  + qez*dpe  - qezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    qeze  = ( qezpe  - qez  )/dpe;
 
elseif str2double(bnd(1))==2
    qeze  = ( qez  - qezme  )/dme;
end
 
if str2double(bnd(2))==1
    qezz  = ( qezpz*dmz  - qez*dmz  + qez*dpz  - qezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    qezz  = ( qezpz  - qez  )/dpz;
 
elseif str2double(bnd(2))==2
    qezz  = ( qez  - qezmz  )/dmz;
end

%% model equations:
section = "Model";
% load model specific equations defined in model.m-file
model

%% 
CC = [];
for icc=1:ncc
    eval(['CC = [CC;eq',cstv{icc},'*diff(eq',cstv{icc},',',cstv{icc},')];'])
end
 
F = []; 
for ix=1:nx
    eval(['F = [F;eq',vars{ix},' + Vzero];']) 
end
 
dF = [];
for ix=1:nx
    eval(['dF = [dF diff(F,',vars{ix},')];'])
end
 
X_ = [];
for ix=1:nx_
    eval(['X_ = [X_;',vars_{ix},' + Vzero];'])
end
 
L = [];
for il=1:nl
    eval(['L = [L;',last{il},' + Vzero];'])
end
 
dFvars = symvar(dF);
Fvars  = symvar(F);
Xvars  = symvar(X_);
Lvars  = symvar(L);
CCvars = symvar(CC);
 
cd('tmp')
eval(['matlabFunction(dF,''File'',''dFfun',bnd,'_'',''Optimize'',false,''Vars'',dFvars);'])
eval(['matlabFunction(F,''File'',''Ffun',bnd,'_'',''Optimize'',false,''Vars'',Fvars);'])
eval(['matlabFunction(X_,''File'',''Xfun',bnd,'_'',''Optimize'',false,''Vars'',Xvars);'])
eval(['matlabFunction(L,''File'',''Lfun',bnd,'_'',''Optimize'',false,''Vars'',Lvars);'])
if ~isempty(CC)
    eval(['matlabFunction(CC,''File'',''CCfun',bnd,'_'',''Optimize'',false,''Vars'',CCvars);'])
end
cd('..')
 

