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
 

if sum(ismember('q', vars)) == 1
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
else
end 

if sum(ismember('qa', vars)) == 1
%% System of Equations for qa
if str2double(bnd(1))==1 
    qae   = ( qape*dme   - qa*dme   + qa*dpe   - qame*dpe   )/(2*dpe*dme);
    qaee  = ( qaepe*dme  - qae*dme  + qae*dpe  - qaeme*dpe  )/(2*dpe*dme);
    qaeee = ( qaeepe*dme - qaee*dme + qaee*dpe - qaeeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    qae   = qaepe - qaeepe*dpe;
    qaee  = qaeepe - qaeeepe*dpe;
    qaeee = ( qaeepe - qaee )/dpe;
 
elseif str2double(bnd(1))==2
    qae   = qaeme + qaeeme*dme;
    qaee  = qaeeme + qaeeeme*dme;
    qaeee = ( qaee - qaeeme )/dme;

end
 
if str2double(bnd(2))==1       
    qaz   = ( qapz*dmz   - qa*dmz   + qa*dpz   - qamz*dpz  )/(2*dpz*dmz);
    qazz  = ( qazpz*dmz  - qaz*dmz  + qaz*dpz  - qazmz*dpz )/(2*dpz*dmz);
    qazzz = ( qazzpz*dmz - qazz*dmz + qazz*dpz - qazzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    qaz   = qazpz - qazzpz*dpz;
    qazz  = qazzpz - qazzzpz*dpz;
    qazzz = ( qazzpz - qazz )/dpz;
    
elseif str2double(bnd(2))==2
    qaz   = qazmz  + qazzmz*dmz;
    qazz  = qazzmz + qazzzmz*dmz;
    qazzz = ( qazz - qazzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    qaze_ = ( qazpe*dme - qaz*dme + qaz*dpe - qazme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    qaze_ = qaezpe - qaezepe*dpe;
    
elseif str2double(bnd(1))==2
    qaze_ = qaezme + qaezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    qaez_ = ( qaepz*dmz - qae*dmz + qae*dpz - qaemz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    qaez_  = qaezpz - qaezzpz*dpz;
 
elseif str2double(bnd(2))==2
    qaez_  = qaezmz + qaezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    qae   = 0;
    qaee  = 0;
    qaeee = 0;
    
    qaz   = 0;
    qazz  = 0;
    qazzz = 0;
    
    qaez  = 0;
    qaeze = 0;
    qaezz = 0;
 
else
    qaez = 1/2*qaez_ + 1/2*qaze_;
end
 
if str2double(bnd(1))==1
    qaeze  = ( qaezpe*dme  - qaez*dme  + qaez*dpe  - qaezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    qaeze  = ( qaezpe  - qaez  )/dpe;
 
elseif str2double(bnd(1))==2
    qaeze  = ( qaez  - qaezme  )/dme;
end
 
if str2double(bnd(2))==1
    qaezz  = ( qaezpz*dmz  - qaez*dmz  + qaez*dpz  - qaezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    qaezz  = ( qaezpz  - qaez  )/dpz;
 
elseif str2double(bnd(2))==2
    qaezz  = ( qaez  - qaezmz  )/dmz;
end

else
end

if sum(ismember('qb', vars)) == 1
%% System of Equations for qb
if str2double(bnd(1))==1 
    qbe   = ( qbpe*dme   - qb*dme   + qb*dpe   - qbme*dpe   )/(2*dpe*dme);
    qbee  = ( qbepe*dme  - qbe*dme  + qbe*dpe  - qbeme*dpe  )/(2*dpe*dme);
    qbeee = ( qbeepe*dme - qbee*dme + qbee*dpe - qbeeme*dpe )/(2*dpe*dme);

elseif str2double(bnd(1))==0
    qbe   = qbepe - qbeepe*dpe;
    qbee  = qbeepe - qbeeepe*dpe;
    qbeee = ( qbeepe - qbee )/dpe;
 
elseif str2double(bnd(1))==2
    qbe   = qbeme + qbeeme*dme;
    qbee  = qbeeme + qbeeeme*dme;
    qbeee = ( qbee - qbeeme )/dme;

end
 
if str2double(bnd(2))==1       
    qbz   = ( qbpz*dmz   - qb*dmz   + qb*dpz   - qbmz*dpz  )/(2*dpz*dmz);
    qbzz  = ( qbzpz*dmz  - qbz*dmz  + qbz*dpz  - qbzmz*dpz )/(2*dpz*dmz);
    qbzzz = ( qbzzpz*dmz - qbzz*dmz + qbzz*dpz - qbzzmz*dpz )/(2*dpz*dmz);
     
elseif str2double(bnd(2))==0
    qbz   = qbzpz - qbzzpz*dpz;
    qbzz  = qbzzpz - qbzzzpz*dpz;
    qbzzz = ( qbzzpz - qbzz )/dpz;
    
elseif str2double(bnd(2))==2
    qbz   = qbzmz  + qbzzmz*dmz;
    qbzz  = qbzzmz + qbzzzmz*dmz;
    qbzzz = ( qbzz - qbzzmz )/dmz;
 
end
 
if str2double(bnd(1))==1    
    qbze_ = ( qbzpe*dme - qbz*dme + qbz*dpe - qbzme*dpe)/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    qbze_ = qbezpe - qbezepe*dpe;
    
elseif str2double(bnd(1))==2
    qbze_ = qbezme + qbezeme*dme;
    
end
 
if str2double(bnd(2))==1        
    qbez_ = ( qbepz*dmz - qbe*dmz + qbe*dpz - qbemz*dpz)/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    qbez_  = qbezpz - qbezzpz*dpz;
 
elseif str2double(bnd(2))==2
    qbez_  = qbezmz + qbezzmz*dmz;
 
end
 
if and(str2double(bnd(1))==3,str2double(bnd(2))==3)
    qbe   = 0;
    qbee  = 0;
    qbeee = 0;
    
    qbz   = 0;
    qbzz  = 0;
    qbzzz = 0;
    
    qbez  = 0;
    qbeze = 0;
    qbezz = 0;
 
else
    qbez = 1/2*qbez_ + 1/2*qbze_;
end
 
if str2double(bnd(1))==1
    qbeze  = ( qbezpe*dme  - qbez*dme  + qbez*dpe  - qbezme*dpe  )/(2*dpe*dme);
 
elseif str2double(bnd(1))==0
    qbeze  = ( qbezpe  - qbez  )/dpe;
 
elseif str2double(bnd(1))==2
    qbeze  = ( qbez  - qbezme  )/dme;
end
 
if str2double(bnd(2))==1
    qbezz  = ( qbezpz*dmz  - qbez*dmz  + qbez*dpz  - qbezmz*dpz  )/(2*dpz*dmz);
 
elseif str2double(bnd(2))==0
    qbezz  = ( qbezpz  - qbez  )/dpz;
 
elseif str2double(bnd(2))==2
    qbezz  = ( qbez  - qbezmz  )/dmz;
end

else
end

%% model equations:
mod_model

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
 
cd(par.tmpFolder)
eval(['matlabFunction(dF,''File'',''dFfun',bnd,'_'',''Optimize'',false,''Vars'',dFvars);'])
eval(['matlabFunction(F,''File'',''Ffun',bnd,'_'',''Optimize'',false,''Vars'',Fvars);'])
eval(['matlabFunction(X_,''File'',''Xfun',bnd,'_'',''Optimize'',false,''Vars'',Xvars);'])
eval(['matlabFunction(L,''File'',''Lfun',bnd,'_'',''Optimize'',false,''Vars'',Lvars);'])
if ~isempty(CC)
    eval(['matlabFunction(CC,''File'',''CCfun',bnd,'_'',''Optimize'',false,''Vars'',CCvars);'])
end
cd('..')
 

