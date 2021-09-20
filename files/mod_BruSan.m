%% PARAMETERS
param = {'gammai','gammah','ai','ah','rhoi','rhoh','sigz','sigbar','deltai','deltah','kappa_p','kappa_z','zetai','zetah'};

par.gammai = 2;
par.gammah = 3;

par.ai     = 0.10;
par.ah     = 0.10;

par.rhoi   = 0.04;
par.rhoh   = 0.04;

par.sigz   = 0.10;
par.sigbar = 0.10;

par.deltai = 0.04;
par.deltah = 0.04;

par.kappa_p = 2;
par.kappa_z = 1;

par.zetai = 1.15;
par.zetah = 1.15;

%% VARIABLES
state  = {'e','z'};
dvec   = {'dme','dmz','dpe','dpz','dae','daz'};
value  = {'vi','vei','vzi','vh','veh','vzh'};

vars   = {'q','psii','mue','sigqk','sigqs'};
prices = {'q'};

vars_  = {'wi','wh','ci','ch',...
          'sigma','siges','sigek',...
          'muri','murh',...
          'muk','muq','r',...
          'muni','munh','signis','signhs','signik','signhk',...
          'sigxis','sigxhs','sigxik','sigxhk','iotai','iotah','muz',...
          'sigV11','sigV22','sigV12','muV1','muV2','uVi','uVh','rVi','rVh'}; 
      
vplot_ = {'vi','vei','vzi','vh','veh','vzh','qz','muz'};

latex_ = {'wi','wh','ci','ch',...
          'sigma','siges','sigek',...
          'muri','murh',...
          'muk','muq','r',...
          'muni','munh','signis','signhs','signik','signhk',...
          '$\sigma^{\xi}$','$\sigma^{\xi}$','$sigxik$','$sigxhk$','$\iota^i$','$\iota^h$','$\mu^z$',...
          'sigV11','sigV22','sigV12','muV1','muV2','uVi','uVh','rVi','rVh'}; 

latex  = {'$q$','$\psi$','$mue$','$sigqk$','$sigqs$'};

latexv = {'$v^i$','$v_\eta^i$','$v_z^i$','$v^h$','$v_\eta^h$','$v_z^h$'};

latexs = {'$\eta$', '$\sigma$'};

%% GUESS
X0 = NaN([par.nx par.dim]);

X0(par.X.iq,:,:) = 1;
X0(par.X.ipsii,:,:) = 0.95;
X0(par.X.isigqs,:,:) = 0;
X0(par.X.isigqk,:,:) = 0;
X0(par.X.imue,:,:) = 0;

vi = par.rhoi*ones(par.dim);
vh = par.rhoh*ones(par.dim);

%% MODEL
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

% volatilities derived from model set up
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
eqmue = ( muni - muk - muq - sigma*sigqk +(sigqk+sigma)^2 + sigqs^2 - wi*sigqs^2 - wi*(sigqk+sigma)^2 ) - mue; % Equation (31)

eqq   = (ci*e + ch*(1-e))*q - psii*(ai-iotai) - (1-psii)*(ah-iotah); % Equation (32)

eqpsii = muri - murh... 
        + gammah*wh*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
        - gammai*wi*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
        + sigqs*sigxis + (sigqk+sigma)*sigxik ...
        - sigqs*sigxhs - (sigqk+sigma)*sigxhk; % Equation (33)
       
eqsigqs = (sigz*z*qz + siges*e*qe) - sigqs*q; % Equation (34)
eqsigqk = sigek*e*qe - sigqk*q; % Equation (35)

% hjb variables
muV1  = mue*e;
muV2  = muz*z;

sigV11 = (siges*e)^2 + (sigek*e)^2;
sigV22 = (sigz*z)^2;
sigV12 =  siges*e*sigz*z;
 
uVi   = 0;
uVh   = 0;

rVi   = - (1-gammai)*( 1/(1-1/zetai)*(ci-rhoi) ...
     + r - ci ...
     + gammai/2*( wi*(sigqs)^2 + wi*(sigqk+sigma)^2 ) ); %Equation (36)

rVh   = - (1-gammah)*( 1/(1-1/zetah)*(ch-rhoh) ...
     + r - ch ...
     + gammah/2.*( wh*(sigqs)^2 + wh*(sigqk+sigma)^2 ) ); %Equation (36)

%% CONSTRAINTS
cstv   = {};
csts   = {};
cstn   = {};

% impose value at boundaries if necessary
par.nbv = zeros(3,3);
bndv = cell(3,3,2);
bndn = cell(3,3,2);
par.base = 3;
