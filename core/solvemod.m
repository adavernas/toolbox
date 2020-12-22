function solvemod(name,varargin)
% This file computes the solution method of d'Avernas and Vandeweyer (2018)
% written by Adrien d'Avernas
% 2019 d'Avernas and Vandeweyer all rights reserved

set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:interp2:NaNstrip') 

if any(strcmp(varargin,'parrallel'))
    default_parrallel = 'off';
elseif ~any(strcmp(varargin,'dimensions'))
    default_parrallel = 'on';
elseif any(strcmp(varargin,'dimensions'))
    ii = find(strcmp(varargin,'dimensions'));
    if strcmp(varargin{ii+1},'2D')
        default_parrallel = 'on';
    elseif strcmp(varargin{ii+1},'1D')
        default_parrallel = 'off';
    end
end

p = inputParser;
addRequired(p,'name',@ischar);
addParameter(p,'write','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'guess','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'dimensions','2D',@(x) any(validatestring(x,{'1D','2D'})));
addParameter(p,'innerplot','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'outerplot','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'allplot','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'debug','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'parrallel',default_parrallel,@(x) any(validatestring(x,{'on','off'})));      
addParameter(p,'loop','off',@(x) any(validatestring(x,{'on','off'})));      
addParameter(p,'search','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'savegraph','off',@(x) any(validatestring(x,{'on','off'})));
addParameter(p,'method','first_order',@(x) any(validatestring(x,{'first_order','forward'})));
addParameter(p,'dispT',10,@isinteger);
addParameter(p,'HJBupdate','on',@(x) any(validatestring(x,{'on','off'})));

parse(p,name,varargin{:});

par = p.Results;
if strcmp(par.debug,'on')
    dbstop if error
    dbstop if warning
elseif strcmp(par.debug,'off')
    dbclear if error
    dbclear if warning
end
    
%% Write Functions
writemod
par_solver
mod_parameters
mod_variables
mod_constraints

par.dim  = [par.n1 par.n2];
par.dim_ = par.dim-2*par.nextr;
par.ndim = size(par.dim,2);
par.N    = prod(par.dim);
par.N_   = prod(par.dim_);

setupfunctions

%% Initial Guess for X
par.solFolder = ['sol_',par.name];
if strcmp(par.guess,'on')
    cd(par.solFolder)
    guess = load('guess.mat');
    cd ..
    
    X0 = NaN([par.nx par.dim]);
    for ix=1:par.nx
        ixg(1) = find(strcmp(vars(ix),guess.vars));
        if ~isempty(ixg)
            X0tmp = interp2(guess.x2,guess.x1,squeeze(guess.X0(ixg(1),:,:)),x2,x1,'spline');
            X0(ix,:,:) = X0tmp;
        end
    end
    
    X_ = zeros([par.nx_ par.dim]);
    for ix_=1:par.nx_
        ixg = find(strcmp(vars_(ix_),guess.vars_));
        if ~isempty(ixg)
            X_tmp = interp2(guess.x2,guess.x1,squeeze(guess.X_(ixg(1),:,:)),x2,x1,'spline');
            X_(ix_,:,:) = X_tmp;
        end
    end
    
    V0 = NaN([par.nv par.dim]);
    for iv=1:par.nv
        ixg(1) = find(strcmp(value(iv),guess.value));
        if ~isempty(ixg)
            V0tmp = interp2(guess.x2,guess.x1,squeeze(guess.V0(ixg(1),:,:)),x2,x1,'spline');
            V0(iv,:,:) = V0tmp;
        end
    end
    
    L = zeros([par.nl par.dim]);
    for il=1:par.nl
        ixg = find(strcmp(last(il),guess.last));
        if ~isempty(ixg)
            Ltmp = interp2(guess.x2,guess.x1,squeeze(guess.L(ixg(1),:,:)),x2,x1,'spline');
            L(il,:,:) = Ltmp;
        end
    end
else
    % call initial guess specified in model.m-file   
    mod_guess

    % compute initial guesses for the partial derivatives of the value functions
    vei = finitediff2D(vi,D(1,:,:),x1,x2,1); %#ok<NODEF>
    vzi = finitediff2D(vi,D(2,:,:),x1,x2,2);
    
    veh = finitediff2D(vh,D(1,:,:),x1,x2,1); %#ok<NODEF>
    vzh = finitediff2D(vh,D(2,:,:),x1,x2,2);
    
    V0 = NaN([par.nv par.dim]);
    V0(par.V.ivi,:,:)  = vi;
    V0(par.V.ivei,:,:) = vei;
    V0(par.V.ivzi,:,:) = vzi;
    V0(par.V.ivh,:,:)  = vh;
    V0(par.V.iveh,:,:) = veh;
    V0(par.V.ivzh,:,:) = vzh;
    
    X_ = zeros([par.nx_ par.dim]);
    
    L = zeros([par.nl par.dim]);
end

%% Iterate over Time for Boundaries
tic

% OUTER LOOP
iterT = 0;
testT = 1;
loop  = 1;
par.minT = par.minT1;
while or(iterT<par.minT,and(testT>=par.tolT,iterT<par.maxT))
iterT = iterT+1;
if round(iterT/par.dispT)==iterT/par.dispT
    msg = ['loop ',num2str(loop),' iteration ',num2str(iterT,'%0.0f'),' with error ',num2str(testT,'%0.6f'),' and dt = ',num2str(par.dt,'%0.2f')];
    disp(msg)
end

iterK = 0;
testK = 1;

while or(iterK<par.minK,and(testK>=par.tolK,iterK<par.maxK))
iterK = iterK+1;
for il=1:par.nl
if  any(strcmp(vars,last{il}))
    eval(['L(par.L.i',last{il},',:,:)  = X0(par.X.i',last{il},',:,:);'])
end
end
for il=1:par.nl
if  and(not(any(strcmp(vars,last{il}))),any(strcmp(vars_,last{il})))
    eval(['L(par.L.i',last{il},',:,:)  = X_(par.X_.i',last{il},',:,:);'])
end
end
for il=1:par.nl
if  any(strcmp(value,last{il}))
    eval(['L(par.L.i',last{il},',:,:)  = V0(par.V.i',last{il},',:,:);'])
end
end

XT0 = X0;
L0 = L;

% INNER LOOP: the following part of code solves the set of non-linear equations for the inner loop
if strcmp(par.parrallel,'off')
    innerloop
elseif strcmp(par.parrallel,'on')
    innerloop_par
end

% once inner loop is done, the function HJB computes the coefficients in bellman equations using the values
% from inner loop

if loop==1
    X0 = (1-par.damp1)*X0 + par.damp1*XT0;
end
if loop==2
    X0 = (1-par.damp2)*X0 + par.damp2*XT0;
end

if strcmp(par.dimensions,'1D')
    for ik = [1:par.nextr par.nextr+2:par.dim(2)]
        X0(:,:,ik) = X0(:,:,par.nextr+1);
        V0(:,:,ik) = V0(:,:,par.nextr+1);
        L(:,:,ik)  = L(:,:,par.nextr+1);
        L0(:,:,ik) = L0(:,:,par.nextr+1);
        P1(:,:,ik) = P1(:,:,par.nextr+1);
        P2(:,:,ik) = P2(:,:,par.nextr+1);
        M1(:,:,ik) = M1(:,:,par.nextr+1);
        M2(:,:,ik) = M2(:,:,par.nextr+1);
    end
end

xxfun

if strcmp(par.innerplot,'on')
    par.savegraph_ = 'off';
    plotgraphs2D
    drawnow
end

testK = sum(sum(sum( (XT0 - X0).^2 )));

end

if loop==1
    par.damp1 = par.damp1_*par.damp1;
    par.damp1 = max(par.damp1,par.dampm1);
end
if loop==2
    par.damp2 = par.damp2_*par.damp2;
    par.damp2 = max(par.damp2,par.dampm2);
end

% stencil decomposition happens here to solve for the PDEs; the coefficients for partial derivaties are computed
% from values in inner loop

vvi = HJBmap2D(DD,VV(par.V.ivi,:),XX_(par.X_.irVi,:),XX_(par.X_.iuVi,:),...
    [XX_(par.X_.imuV1,:);XX_(par.X_.imuV2,:)],...
    [XX_(par.X_.isigV11,:);XX_(par.X_.isigV22,:);XX_(par.X_.isigV12,:)],par);
vvh = HJBmap2D(DD,VV(par.V.ivh,:),XX_(par.X_.irVh,:),XX_(par.X_.iuVh,:),...
    [XX_(par.X_.imuV1,:);XX_(par.X_.imuV2,:)],...
    [XX_(par.X_.isigV11,:);XX_(par.X_.isigV22,:);XX_(par.X_.isigV12,:)],par);

vi  = squeeze(reshape(vvi,[1 par.dim]));  
vh  = squeeze(reshape(vvh,[1 par.dim]));    

[vei,~,~] = finitediff2D(vi,D(1,:,:),x1,x2,1);
[vzi,~,~] = finitediff2D(vi,D(2,:,:),x1,x2,2);
[veh,~,~] = finitediff2D(vh,D(1,:,:),x1,x2,1);
[vzh,~,~] = finitediff2D(vh,D(2,:,:),x1,x2,2);

% compute the value function numerical error testT; this is compared against the tolerance level in outerloop
testT = max(max( abs( ( vi-squeeze(V0(par.V.ivi,:,:)) )./squeeze(V0(par.V.ivi,:,:)) )));
  
if strcmp(par.loop,'on')
if loop==1
if and(testT<par.tolT,iterT>=par.minT)
    loop  = 2;
    testT = 1;
    iterT = 0;
    par.dt   = par.dt2;
    par.minT = par.minT2;
end
end
end

if strcmp(par.HJBupdate,'on')
    V0(par.V.ivi,:,:)  = vi;
    V0(par.V.ivei,:,:) = vei;
    V0(par.V.ivzi,:,:) = vzi;
    V0(par.V.ivh,:,:)  = vh;
    V0(par.V.iveh,:,:) = veh;
    V0(par.V.ivzh,:,:) = vzh;
end 
        
if round(iterT/par.dispT)==iterT/par.dispT
if strcmp(par.outerplot,'on')
    par.savegraph_ = 'off';
    plotgraphs2D    
    drawnow
end
end

end

toc

if ~exist(par.solFolder, 'dir')
    mkdir(par.solFolder)
end
cd(par.solFolder)
tmp = 'sol';
for i=1:length(param)
    eval(['tmp = [tmp,''_'',''',param{i},''',''_'',num2str(round(par.',param{i},'*100))];'])
end
save(tmp,'X0','X_','V0','L','vars','vars_','value','last','x1','x2','vec1','vec2','par')
save('guess','X0','X_','V0','L','vars','vars_','value','last','x1','x2','vec1','vec2','par')
cd ..

par.savegraph_ = par.savegraph;
plotgraphs2D  
