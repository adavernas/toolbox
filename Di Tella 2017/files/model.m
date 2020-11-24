%% Define model specific parameters:
if (section == "Parameters")
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
        

%% Define model specific variables:
elseif (section == "Variables")

vars   = {'q','mue','sigw', 'sigx','theta'};

% all secondary variables are put in vars_ MATLAB variable; these secondary variables are computed from
% the endogenous variables in vars; the mapping is given in the paper and
% in writefun.m
vars2  = {'ci','ch',...
          'sigz_',...
          'muq','r',...
          'muni','munh','sign',...
          'sigxii','sigxih','g','muz', ...
          'pii', 'iotai', 'v', 'x', 'k', 'sigq'};
      
% these are the latex symbols used for graphs
latex  = {'$q$','$\mu^{\eta}$','$\sigma^w$','$\sigma^x$','$\theta$'};


latex2 ={'ci','ch',...
          'sigma','sigz_',...
          'muq','r',...
          'muni','munh','sign',...
          'sigxii','sigxih','g','muz', ...
          'pii', 'iotai', 'v', 'x', 'k', 'sigq'};

latexv = {'$v^i$','$v_\eta^i$','$v_z^i$','$v^h$','$v_\eta^h$','$v_z^h$'};

%% Initial guess for variables:
elseif(section == "Guess") && strcmp(par.guess,'off')
    %if no initial exogenous guess file is provided, initial values should
    %be given here. These are also model specific and have to be adapted to
    %the model variables. X0 contains the main model variables, while secondary variables are in X_ 
X0 = NaN([par.nx par.dim]);
    
    X0(par.X.iq,:,:) = 1;   
    X0(par.X.isigw,:,:) = eps;
    X0(par.X.isigx,:,:) = eps;
    X0(par.X.itheta,:,:) = eps;
    
    X0(par.X.imue,:,:) = eps;
    
    %X0(par.X.ipii,:,:) = eps;
    
    X_(:,:,:)           = eps;
    %X_(par.X_.iqe,:,:) = eps;
    %X_(par.X_.iqz,:,:) = eps;
    
    %X_(par.X_.iqee,:,:) = eps;
    %X_(par.X_.iqzz,:,:) = eps;
        
    % Initial Guess for V
    vi_ = par.rhoi;
    vh_ = par.rhoh;
    
    
elseif (section == "Model")
%% Model specific set up:
%state variables in this model are eta (e) and sigma (z):
v = z; 
x = e;
%foc wrt to consumption:
ci    = rhoi^(1/psii)*vi^((psii - 1)/psii);  
ch    = rhoi^(1/psii)*vh^((psii - 1)/psii); 

k = 1/(q*x);

%investment maximization:
g  = 1/(2*A)*(q - B) - deltai; %A and B calibrated such that GDP growth is equal to 2%
% per setup of the model:
iotai  = A*(g + deltai)^2 + B*(g + deltai) ;

muz   = kappa_z*(vbar - v);
sigz_ = sigz*sqrt(v)        ;

    
% Volatilities 
% value function is a function of the two state variables sigma (z) and eta
% (e)
sigxii = vei/vi*sigx*e + vzi/vi*sigz_;
sigxih = veh/vh*sigx*e + vzh/vh*sigz_;

sigq  = qe/q*sigx   + qz/q*sigz_*sqrt(z) ;

sign  = phii_*q*k*(sigma + sigq); %sigma + sigq + sigx/e;

pii   = gammai*sigw + (gammai - 1)*sigxih;


muq   =  qe/q*mue*e + qz/q*muz*z ...
         + 1/2*qee/q*( sigx^2 )...
         + 1/2*qzz/q*(sigz_)^2*z ...
         + qez/q*sigx*sigz_*sqrt(z);       

r = (ai - iotai)/q + g + muq + sigma*sigq - ( 1 - phii_)*(sigma + sigq)*pii ...
    - gammai*(sigma + sigq)*(phii_*q*k*(sigma + sigq) + theta) ...
    + ( 1 - gammai)*phii_*(sigma + sigq)*sigxii ...
    - gammai*q*k*(phii_*z)^2;

muni  = r + gammai*(phii_*z)^2/e^2 + pii*sign;
munh  = r + pii*sigw; 




% System to Solve with Newton-Raphson
eqmue = (muni - ci - muq - g - sigq*sigma + (sigma + sigq)^2 - sign*(sigma + sigq))*x - mue; 
    
eqq   = (ci*e + ch*(1-e))*q - (ai-g); 

eqsigw = sign*e + sigw*(1 - e) - sigma - sigq;
eqsigx = (sign + sigma + sigq)*e - sigx ; 
eqtheta = pii + ( 1 - gammai)*sigxii - gammai*phii_*q*k*(sigma + sigq) - gammai*theta;

end


