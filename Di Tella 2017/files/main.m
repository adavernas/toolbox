% This file computes the solution method of d'Avernas and Vandeweyer (2018)
% written by Adrien d'Avernas
% 2019 d'Avernas and Vandeweyer all rights reserved

clearvars
set(0,'DefaultFigureWindowStyle','docked')
cdtmp = cd;
cd ..
addpath([pwd,'/core'])
cd(cdtmp)
addpath([pwd,'/tmp'])

warning('off','MATLAB:interp2:NaNstrip')

dbstop if error
dbstop if warning

par.name = 'davan2018';  % save file name as davan2018
par.dispT = 1;
par.showgraphs1 = 'off'; % show graphs in first loop
par.showgraphs2 = 'off'; % show graphs in second loop
par.write       = 'on'; % if you change mdoel equations in model.m, change par.write to on and then run the code
par.write_con   = 'on';
par.plotall     = 'off';
par.savegraph   = 'off';
par.guess       = 'off';
par.loop        = 'on';

% ALGORITHM OPTIONS
% par.alg = 'GradientDescent'; % algorithm to solve the non-linear set of equations in inner loop
par.alg = 'NewtonRaphson';

% STATIC LOOP
par.minI = 2;    % necessary for the constraints to have time to adjust
par.maxI = 10000;
par.tolI = 1e-6; % make sure that par.tolI is not higher than par.tolK

% OUTER STATIC LOOP: to update the derivatives using past iterations
par.minK = 2;
par.maxK = 100;
par.tolK = 1e-4;

par.damp1 = 0.0;
par.damp1_ = 0.0;
par.dampm1 = 0.0;

par.damp2 = 0.0;
par.damp2_ = 0.0;
par.dampm2 = 0.0;

% HJB LOOP: to update the derivatives using past iterations
par.minT1 = 10;
par.minT2 = 30;

par.maxT = 10000;
par.tolT = 0.01;

% with implicit method, dt can be large; however, it cannot be too big since we use guesses in the static
% loop; if there are convergence issues, reduce dt and try again; dt2 should be small
par.dt  = 1; 
par.dt2 = 1;

% STENCIL
par.maxD = 1000;
par.tolD = 1e-4;

% maxP refers to how far we go in the direction in stencil; refer to the algorithm in paper for details
par.maxP  = 2;
par.nextr = par.maxP; % taille de l'extrapolation

% STATE VARIABLES GRID
vec1a = 0.01;
vec1b = 0.99;
vec2a = 0.01;
vec2b = 0.99;

par.n1 = 24;
par.n2 = 24; 

vec1 = linspace(vec1a,vec1b,par.n1); % linear spacing is generally more stable
vec2 = linspace(vec2a,vec2b,par.n2);

par.dim  = [par.n1 par.n2];
par.dim_ = par.dim-2*par.nextr;
par.ndim = size(par.dim,2);
par.N    = prod(par.dim);
par.N_   = prod(par.dim_);

%% PARAMETERS
section = "Parameters";
%call model specific parameters defined in model.m-file
model
%% DECLARE VARIABLES
section = "Variables";
% state variables are e and z. dme, dpe, dae refers to difference to the left, difference to the right, and
% average differnece for state variable e; similarly we have dmz, dpz; vi, vei, vzi, etc refers to the
% value functions and derivative of value functions.
% Typically, state and dvec are not model specific and are not changed.
% value will also mostly be the same
state  = {'e','z'};
dvec   = {'dme','dmz','dpe','dpz','dae','daz'};
value  = {'vi','vei','vzi','vh','veh','vzh'};


% Call model specific variables defined in model.m-file
model

vars_  = {'qe', 'qee', ...
          'qz', 'qzz', ...
          'qez','qeze','qezz'};
vars_  = [vars_, vars2];      

% these contains all variables that we need from the previous iteration
last   = {'q',...
          'qe', 'qee', 'qeee',...
          'qz', 'qzz', 'qzzz',...
          'qez','qeze','qezz'};


latex_ = {'$q_\eta$','$q_{\overline{\eta}}$','$q_{\eta\eta}$','$q_{\overline{\eta}\overline{\eta}}$',...
          '$q_{{\eta}\overline{\eta}}$','$q_{\eta\overline{\eta}\eta}$','$q_{\eta\overline{\eta}\overline{\eta}}$',...
          '$q_\eta$','$q_{\overline{\eta}}$','$q_{\eta\eta}$','$q_{\overline{\eta}\overline{\eta}}$',...
          '$q_{{\eta}\overline{\eta}}$','$q_{\eta\overline{\eta}\eta}$','$q_{\eta\overline{\eta}\overline{\eta}}$'};
latex_ = [latex_, latex2];      
      
%% CONSTRAINTS
%boundary constraints can be specified here
cstv   = {};
csts   = {};
cstn   = {};

cstv_ = {}; 
cstn_ = {};

%% Replace Always if Condition Met
cstva  = {};
cstsa  = {};
cstna  = {};

cstva_ = {};
cstna_ = {};

%% 
conv = {};
cons = {};
conn = {};

%% Impose value at boundaries if necessary
par.nbv = zeros(3,3);
bndv = cell(3,3,2);
bndn = cell(3,3,2);
par.base = 3;
    
%% Write Functions
setupfunctions

%% Initial Guess for X
section = "Guess";
if strcmp(par.guess,'on')
    guess = load('guess.mat');
    
    X0 = NaN([par.nx par.dim]);
    for ix=1:par.nx
        X0tmp = interp2(guess.x2,guess.x1,squeeze(guess.X0(ix,:,:)),x2,x1,'spline');
        X0(ix,:,:) = X0tmp;
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
        V0tmp = interp2(guess.x2,guess.x1,squeeze(guess.V0(iv,:,:)),x2,x1,'spline');
        V0(iv,:,:) = V0tmp;
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
 % call intitial guess specified in model.m-file   
 model
 
 % Compute initial guesses for the partial derivatives of the value
 % functions
 vi = NaN(par.n1,par.n2);
 vh = NaN(par.n1,par.n2);
     for i1=1:par.n1
        for i2=1:par.n2
            vi(i1,i2) = vi_;
            vh(i1,i2) = vh_;
        end
    end
    vei = finitediff2D(vi,D(1,:,:),x1,x2,1);
    vzi = finitediff2D(vi,D(2,:,:),x1,x2,2);
    
    veh = finitediff2D(vh,D(1,:,:),x1,x2,1);
    vzh = finitediff2D(vh,D(2,:,:),x1,x2,2);
    
    V0 = NaN([par.nv par.dim]);
    V0(par.V.ivi,:,:)  = vi;
    V0(par.V.ivei,:,:) = vei;
    V0(par.V.ivzi,:,:) = vzi;
    V0(par.V.ivh,:,:)  = vh;
    V0(par.V.iveh,:,:) = veh;
    V0(par.V.ivzh,:,:) = vzh;
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
LT0 = L;

% INNER LOOP: the following part of code solves the set of non-linear equations for the inner loop
%% Newton-Raphson Method
itervec{1,1} = 1:par.dim(1);
itervec{1,2} = 1:par.dim(2);

for it=1:size(itervec,1)
for i2=itervec{it,2}
for i1=itervec{it,1}
        
    L0 = LT0;
    
    if i1>1
        M1(:,i1,i2) = L0(:,i1-1,i2); %#ok<SAGROW>
    elseif i1==1
        M1(:,i1,i2) = NaN; %#ok<SAGROW>
    end

    if i1<par.n1
        P1(:,i1,i2) = L0(:,i1+1,i2); %#ok<SAGROW>
    elseif i1==par.n1
        P1(:,i1,i2) = NaN; %#ok<SAGROW>
    end
    
    if i2>1
        M2(:,i1,i2) = L0(:,i1,i2-1); %#ok<SAGROW>
    elseif i2==1
        M2(:,i1,i2) = NaN; %#ok<SAGROW>
    end
    
    if i2<par.n2
        P2(:,i1,i2) = L0(:,i1,i2+1); %#ok<SAGROW>
    elseif i2==par.n1
        P2(:,i1,i2) = NaN; %#ok<SAGROW>
    end
    
    if i1>1
        XX0 = X0(:,i1-1,i2);
    else
        XX0 = X0(:,i1,i2);
    end
    
    iterI = 0;
    testI = 1;

    X1    = XX0;
    Stmp  = S(:,i1,i2);
    Dtmp  = D(:,i1,i2);
    V0tmp = V0(:,i1,i2);
    L0tmp = L0(:,i1,i2);
    M1tmp = M1(:,i1,i2);
    M2tmp = M2(:,i1,i2);
    P1tmp = P1(:,i1,i2);
    P2tmp = P2(:,i1,i2);
    CCtmp = NaN;
    COtmp = NaN; 
    
    [vtmp,~,BC,BC_,icst,icst_] = verfun(Stmp,XX0,NaN(par.ncc,1),cstv,csts,cstn,bndv,bndn,i1,i2,par);
      
    if strcmp(par.loop,'on')
    if loop==1
        vtmp = 0;
    end
    end
        
    ctmp1 = zeros(par.ncc,1);
    while or(iterI<par.minI,and(testI>=par.tolI,iterI<par.maxI))
        iterI = iterI+1;

        process1_
        
        ix = and(not(ismember(1:par.nx,icst(ctmp1==1))),not(ismember(1:par.nx,icst_(ctmp1==1))));
        
        if rcond(dFF(ix,ix))<1e-16
            disp(['rank not full at ',num2str(i1),' ',num2str(i2)])
        end
        
        if or(iterI>10,strcmp(par.alg,'NewtonRaphson'))
            % Precondition Inversion
            XdFF = diag(1./sum(abs(dFF(ix,ix)),1));
            IdFF = XdFF/(dFF(ix,ix)*XdFF);
            
            X1(ix) = XX0(ix) - IdFF*FF(ix);
            
        elseif strcmp(par.alg,'GradientDescent')
            if iterI>1
                nGD = norm(dFF(ix,ix)'*FF(ix)-dFF_(ix,ix)'*FF_(ix));
                if nGD==0
                    gam = 0.1;
                else
                    gam = (XX0(ix)-XX0_(ix))'*(dFF(ix,ix)'*FF(ix)-dFF_(ix,ix)'*FF_(ix))/norm(dFF(ix,ix)'*FF(ix)-dFF_(ix,ix)'*FF_(ix));
                end
            else
                gam = 0.1;
            end
            
            X1(ix) = XX0(ix) - gam*dFF(ix,ix)'*FF(ix);
        end
               
        if iterI>par.minI-1
            [vtmp,ctmp1,BC,BC_,icst,icst_] = verfun(Stmp,X1,NaN(par.ncc,1),cstv,csts,cstn,bndv,bndn,i1,i2,par);
        end
        
        if strcmp(par.loop,'on')
        if loop==1
            vtmp = 0;
        end
        end
      
        X1(icst(ctmp1==1))  = BC(ctmp1==1);
        X1(icst_(ctmp1==1)) = BC_(ctmp1==1);
        
        testI = norm(FF(ix));
        
        dFF_ = dFF;
        FF_  = FF;
        XX0_ = XX0;
        XX0  = real(X1);

    end
    
    process2_
    
    X0(:,i1,i2) = XX0;
end
end
end

% once inner loop is done, the function HJB computes the coefficients in bellman equations using the values
% from inner loop

if loop==1
    X0 = (1-par.damp1)*X0 + par.damp1*XT0;
end
if loop==2
    X0 = (1-par.damp2)*X0 + par.damp2*XT0;
end

xxfun

if strcmp(par.showgraphs1,'on')
    par.savegraph_ = 'off';
    plotgraphs2D
    drawnow
end

testK = sum(sum(sum( (XT0 - X0).^2 )));

end

if strcmp(par.loop,'on')
if loop==1
    par.damp1 = par.damp1_*par.damp1;
    par.damp1 = max(par.damp1,par.dampm1);
end
end

if strcmp(par.loop,'on')
if loop==2
    par.damp2 = par.damp2_*par.damp2;
    par.damp2 = max(par.damp2,par.dampm2);
end
end

[rrvi,uuvi,mmuvi,ssigvi] = HJB(SS,VV,XX,XX_,par,'i'); 
[rrvh,uuvh,mmuvh,ssigvh] = HJB(SS,VV,XX,XX_,par,'h'); 

% stencil decomposition happens here to solve for the PDEs; the coefficients for partial derivaties are computed
% from values in inner loop

vvi = HJBmap2D(DD,VV(par.V.ivi,:),rrvi,uuvi,mmuvi,ssigvi,par);
vvh = HJBmap2D(DD,VV(par.V.ivh,:),rrvh,uuvh,mmuvh,ssigvh,par);

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

V0(par.V.ivi,:,:)  = vi;
V0(par.V.ivei,:,:) = vei;
V0(par.V.ivzi,:,:) = vzi;
V0(par.V.ivh,:,:)  = vh;
V0(par.V.iveh,:,:) = veh;
V0(par.V.ivzh,:,:) = vzh;

if round(iterT/par.dispT)==iterT/par.dispT
if strcmp(par.showgraphs2,'on')
    par.savegraph_ = 'off';
    plotgraphs2D    
    drawnow
end
end

end

tmp = 'sol';
for i=1:length(param)
    eval(['tmp = [tmp,''_'',''',param{i},''',''_'',num2str(round(par.',param{i},'*100))];'])
end


toc

figure(1); clf(1);
par.savegraph_ = par.savegraph;