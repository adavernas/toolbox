% This file specifies the model introduced by Di Tella (2016)
% It is part of the toolbox developed in d'Avernas, Schubert and Vandeweyer (2020)
% The model presented here serves as working example of this toolbox
% written by Valentin Schubert
% 2020 d'Avernas, Schubert and Vandeweyer all rights reserved

%% PARAMETERS
param = {'gammai','ai','ah','rhoi','rhoh','sigz','vbar','deltai','deltah',...
         'kappa_p','kappa_z','zetai','zetah', 'psii','phii_', 'tauu', 'lambda', 'sigma', ...
         'A', 'B'};
par.gammai = 5;


par.sigz   = 0.05;
par.vbar = 0.10;

par.ai     = 0.10;
par.ah     = 0.10;

par.rhoi   = 0.0665;
par.rhoh   = 0.0665;


par.deltai = 0.04;
par.deltah = 0.04;

par.kappa_p = 2;
par.kappa_z = 1;

par.zetai = 1.01;
par.zetah = 1.01;

par.psii  = 1/1.01;

par.phii_  = 0.1 ;
par.tauu   = 0.1 ;
par.lambda = 1  ;

par.A      = 5 ;
par.B      = 2 ;

par.sigma  = 0.0125; %aggregate Brownian, constant in the model
        

%% VARIABLES

state  = {'e','z'};
dvec   = {'dme','dmz','dpe','dpz','dae','daz'};
value  = {'vi','vei','vzi','vh','veh','vzh'};


vars   = {'q','mue','sigw', 'sigx','theta'};

% all secondary variables are put in vars_ MATLAB variable; these secondary variables are computed from
% the endogenous variables in vars; the mapping is given in the paper and
% in writefun.m
vars_  = {'qe', 'qee', ...
          'qz', 'qzz', ...
          'qez','qeze','qezz',...
          'ci','ch',...
          'sigz_',...
          'muq','r',...
          'muni','munh','sign',...
          'sigxii','sigxih','g','muz', ...
          'pii', 'iotai', 'v', 'x', 'k', 'sigq',...
          'sigV11','sigV22','sigV12','muV1','muV2','uVi','uVh','rVi','rVh'};

vplot_ = {'vi'};

last   = {'q',...
          'qe', 'qee', 'qeee',...
          'qz', 'qzz', 'qzzz',...
          'qez','qeze','qezz'};
% these are the latex symbols used for graphs
latex  = {'$q$','$\mu^{\eta}$','$\sigma^w$','$\sigma^x$','$\theta$'};


latex_ ={'$q_\eta$','$q_{\overline{\eta}}$','$q_{\eta\eta}$','$q_{\overline{\eta}\overline{\eta}}$',...
          '$q_{{\eta}\overline{\eta}}$','$q_{\eta\overline{\eta}\eta}$','$q_{\eta\overline{\eta}\overline{\eta}}$',...
          '$q_\eta$','$q_{\overline{\eta}}$','$q_{\eta\eta}$','$q_{\overline{\eta}\overline{\eta}}$',...
          '$q_{{\eta}\overline{\eta}}$','$q_{\eta\overline{\eta}\eta}$','$q_{\eta\overline{\eta}\overline{\eta}}$',...
          'ci','ch',...
          'sigma','sigz_',...
          'muq','r',...
          'muni','munh','sign',...
          'sigxii','sigxih','g','muz', ...
          'pii', 'iotai', 'v', 'x', 'k', 'sigq'};

latexv = {'$v^i$','$v_\eta^i$','$v_z^i$','$v^h$','$v_\eta^h$','$v_z^h$'};




%% GUESS
X0 = NaN([par.nx par.dim]);

X0(par.X.iq,:,:) = 1;   
X0(par.X.isigw,:,:) = eps;
X0(par.X.isigx,:,:) = eps;
X0(par.X.itheta,:,:) = eps;
X0(par.X.imue,:,:) = eps;

X_(par.X_.iqe,:,:) = eps;
X_(par.X_.iqz,:,:) = eps;    
X_(par.X_.iqee,:,:) = eps;
X_(par.X_.iqzz,:,:) = eps;
    

vi = par.rhoi*ones(par.dim);
vh = par.rhoh*ones(par.dim);

%% MODEL
%state variables in this model are eta (e) and sigma (z):
v = z; %Equation (1)
x = e; %Equation (2)
%foc wrt to consumption:
ci    = rhoi^(1/psii)*vi^((psii - 1)/psii);  %Equation (5)
ch    = rhoi^(1/psii)*vh^((psii - 1)/psii);  %Equation (6)

k = 1/(q*x); %Equation (7)

%investment maximization:
g  = 1/(2*A)*(q - B) - deltai;  %Equation (8) //A and B calibrated such that GDP growth is equal to 2%
% per setup of the model:
iotai  = A*(g + deltai)^2 + B*(g + deltai) ;  %Equation (9)

muz   = kappa_z*(vbar - v);  %Equation (10)
sigz_ = sigz*sqrt(v)        ;  %Equation (11)

    
% Volatilities 
% value function is a function of the two state variables sigma (z) and eta
% (e)
sigxii = vei/vi*sigx*e + vzi/vi*sigz_;  %Equation (12)
sigxih = veh/vh*sigx*e + vzh/vh*sigz_;  %Equation (13)

sigq  = qe/q*sigx   + qz/q*sigz_*sqrt(z) ;  %Equation (14)

sign  = phii_*q*k*(sigma + sigq); %Equation (16)

pii   = gammai*sigw + (gammai - 1)*sigxih; %Equation (17)


muq   =  qe/q*mue*e + qz/q*muz*z ...  %Equation (15)
         + 1/2*qee/q*( sigx^2 )...
         + 1/2*qzz/q*(sigz_)^2*z ...
         + qez/q*sigx*sigz_*sqrt(z);       

r = (ai - iotai)/q + g + muq + sigma*sigq - ( 1 - phii_)*(sigma + sigq)*pii ...
    - gammai*(sigma + sigq)*(phii_*q*k*(sigma + sigq) + theta) ...
    + ( 1 - gammai)*phii_*(sigma + sigq)*sigxii ...
    - gammai*q*k*(phii_*z)^2; %Equation (18)

muni  = r + gammai*(phii_*z)^2/e^2 + pii*sign; %Equation (19)
munh  = r + pii*sigw; %Equation (20)




% System to Solve with Newton-Raphson
eqmue = (muni - ci - muq - g - sigq*sigma + (sigma + sigq)^2 - sign*(sigma + sigq))*x - mue; %Equation (21)
    
eqq   = (ci*e + ch*(1-e))*q - (ai-g); %Equation (22)

eqsigw = sign*e + sigw*(1 - e) - sigma - sigq; %Equation (23)
eqsigx = (sign + sigma + sigq)*e - sigx ; %Equation (24)
eqtheta = pii + ( 1 - gammai)*sigxii - gammai*phii_*q*k*(sigma + sigq) - gammai*theta; %Equation (25)


% hjb variables
muV1  = mue*e;
muV2  = muz*z;

sigV11 = (sigx*e)^2;
sigV22 = (sigz_*z)^2;
sigV12 =  sigx*e*sigz_*z;
 
uVi   = 0;
uVh   = 0;

rVi  =  - ( 1/(1-psii)*(ci - rhoi) ...
     + r - ci  - theta^2 ...
     + gammai/2*( (sign)^2 )) ;
 
rVh  =  - ( 1/(1-psii)*(ch - rhoh) ...
     + r - ch...
     + gammai/2*( (sigw)^2 )) ; 

%% CONSTRAINTS
cstv   = {};
csts   = {};
cstn   = {};

% impose value at boundaries if necessary
par.nbv = zeros(3,3);
bndv = cell(3,3,2);
bndn = cell(3,3,2);
par.base = 3;
