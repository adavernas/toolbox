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
 
%% System of Equations for derivatives of prices
mod_differences

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
 

