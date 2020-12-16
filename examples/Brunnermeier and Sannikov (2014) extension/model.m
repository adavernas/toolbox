% This file specifies an extension of the model introduced by Brunnermeier and Sannikov (2014)
% It is part of the toolbox developed in d'Avernas, Schubert and Vandeweyer (2020)
% The model presented here serves as working example of this toolbox
% written by Valentin Schubert
% 2020 d'Avernas, Schubert and Vandeweyer all rights reserved

%% Define model specific parameters:
if (section == "Parameters")
param = {'gammai','gammah','ai','ah','rhoi','rhoh','sigz','sigbar','deltai','deltah','kappa_p','kappa_z','zetai','zetah'};

par.gammai = 2;
par.gammah = 3;

par.sigz   = 0.10;
par.sigbar = 0.10;

par.ai     = 0.10;
par.ah     = 0.10;

par.rhoi   = 0.04;
par.rhoh   = 0.04;

par.deltai = 0.04;
par.deltah = 0.04;

par.kappa_p = 2;
par.kappa_z = 1;

par.zetai = 1.01;
par.zetah = 1.01;
%% Define model specific variables:
elseif (section == "Variables")

vars   = {'q','psii','mue','sigqk','sigqs'};

% all secondary variables are put in vars_ MATLAB variable; these secondary variables are computed from
% the endogenous variables in vars; the mapping is given in the paper and
% in writefun.m
vars2  = {'wi','wh','ci','ch',...
          'sigma','siges','sigek',...
          'muri','murh',...
          'muk','muq','r',...
          'muni','munh','signis','signhs','signik','signhk',...
          'sigxis','sigxhs','sigxik','sigxhk','iotai','iotah','muz'};
      
% these are the latex symbols used for graphs
latex  = {'$q$','$\psi$','$mue$','$sigqk$','$sigqs$'};


latex2 = {'phii','phih',...
          'wi','wh','ci','ch',...
          'sigma','siges','sigek',...
          'muri','murh',...
          'muk','muq','r',...
          'muni','munh','signis','signhs','signik','signhk',...
          '$\sigma^{\xi}$','$\sigma^{\xi}$','$sigxik$','$sigxhk$','$\iota^i$','$\iota^h$','$\mu^z$'};

latexv = {'$v^i$','$v_\eta^i$','$v_z^i$','$v^h$','$v_\eta^h$','$v_z^h$'};

%% Initial guess for variables:
elseif(section == "Guess") && strcmp(par.guess,'off')
    %if no initial exogenous guess file is provided, initial values should
    %be given here. These are also model specific and have to be adapted to
    %the model variables. X0 contains the main model variables, while secondary variables are in X_ 
    X0 = NaN([par.nx par.dim]);
    
    X0(par.X.iq,:,:) = 1;   
    X0(par.X.ipsii,:,:) = 0.95;
    X0(par.X.isigqs,:,:) = eps;
    X0(par.X.isigqk,:,:) = eps;
    
    X0(par.X.imue,:,:) = eps;
    
    
    X_(par.X_.iqe,:,:) = eps;
    X_(par.X_.iqz,:,:) = eps;
    
    X_(par.X_.iqee,:,:) = eps;
    X_(par.X_.iqzz,:,:) = eps;
        
    % Initial Guess for V
    vi_ = par.rhoi;
    vh_ = par.rhoh;
    
elseif (section == "Model")
%% Model specific set up:
%state variables in this model are eta (e) and sigma (z):
eta   = e; % Equation (1)
sigma = z; % Equation (2)

wi = psii/e; % Equation (5)
wh = (1-psii)/(1-e); % Equation (6)

%foc wrt to consumption:
ci    = vi^((1-zetai)/(1-gammai)); % Equation (7)
ch    = vh^((1-zetah)/(1-gammah)); % Equation (8)

%investment maximization:
iotai  = (q-1)/kappa_p; % Equation (9)
iotah  = (q-1)/kappa_p; % Equation (10)

% per setup of the model:
phii  = log(1+kappa_p*iotai)/kappa_p - deltai; % Equation (11)
phih  = log(1+kappa_p*iotah)/kappa_p - deltah; % Equation (12)

muz = kappa_z*(sigma - sigbar); % Equation (13)

% from aggregation:
muk   = psii*phii + (1-psii)*phih; % Equation (14)
    
% Volatilities 
% derived from model set up
signis = wi*sigqs; % Equation (15)
signhs = wh*sigqs; % Equation (16)

signik = wi*(sigqk+sigma); % Equation (17)
signhk = wh*(sigqk+sigma); % Equation (18)

siges  = signis -  sigqs       ; % Equation (19)
sigek  = signik - (sigqk+sigma); % Equation (20)

% value function is a function of the two state variables sigma (z) and eta
% (e)
sigxik = vei/vi*sigek*e; % Equation (21)
sigxhk = veh/vh*sigek*e; % Equation (22)

sigxis = vei/vi*siges*e + vzi/vi*sigz*z; % Equation (23)
sigxhs = veh/vh*siges*e + vzh/vh*sigz*z; % Equation (24)

muq   =  qe/q*mue*e + qz/q*muz*z ...
         + 1/2*qee/q*( (siges*e)^2 + (sigek*e)^2 )...
         + 1/2*qzz/q*(sigz*z)^2 ...
         + qez*siges*e*sigz*z;       % Equation (25)

muri  = (ai-iotai)/q + phii + muq + sigma*sigqk; % Equation (26)
murh  = (ah-iotah)/q + phih + muq + sigma*sigqk; % Equation (27)

r = muri ...
    - gammai*wi*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
    + sigqs*sigxis + (sigqk+sigma)*sigxik; % Equation (28)

muni  = r + wi*(muri-r) - ci; % Equation (29)
munh  = r + wh*(murh-r) - ch; % Equation (30)


% System to Solve with Newton-Raphson
eqmue = ( muni - muk - muq - sigma*sigqk + 1/2*(sigqk+sigma)^2 + 1/2*sigqs^2 - wi*sigqs^2 - wi*(sigqk+sigma)^2 ) - mue; % Equation (31)

eqq   = (ci*e + ch*(1-e))*q - psii*(ai-iotai) - (1-psii)*(ah-iotah); % Equation (32)

eqpsii = muri - murh... 
        + gammah*wh*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
        - gammai*wi*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
        + sigqs*sigxis + (sigqk+sigma)*sigxik ...
        - sigqs*sigxhs - (sigqk+sigma)*sigxhk; % Equation (33)
       
eqsigqs = (sigz*z*qz + siges*e*qe) - sigqs*q; % Equation (34)
eqsigqk = sigek*e*qe - sigqk*q; % Equation (35)

end


