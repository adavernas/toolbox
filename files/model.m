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
vars_  = {'qe', 'qee', ...
          'qz', 'qzz', ...
          'qez','qeze','qezz',...
          'wi','wh','ci','ch',...
          'sigma','siges','sigek',...
          'muri','murh',...
          'muk','muq','r',...
          'muni','munh','signis','signhs','signik','signhk',...
          'sigxis','sigxhs','sigxik','sigxhk','iotai','iotah','muz'};

% these contains all variables that we need from the previous iteration
last   = {'q',...
          'qe', 'qee', 'qeee',...
          'qz', 'qzz', 'qzzz',...
          'qez','qeze','qezz'};

% these are the latex symbols used for graphs
latex  = {'$q$','$\psi$','$mue$','$sigqk$','$sigqs$'};

latex_ = {'$q_\eta$','$q_{\overline{\eta}}$','$q_{\eta\eta}$','$q_{\overline{\eta}\overline{\eta}}$',...
          '$q_{{\eta}\overline{\eta}}$','$q_{\eta\overline{\eta}\eta}$','$q_{\eta\overline{\eta}\overline{\eta}}$',...
          '$q_\eta$','$q_{\overline{\eta}}$','$q_{\eta\eta}$','$q_{\overline{\eta}\overline{\eta}}$',...
          '$q_{{\eta}\overline{\eta}}$','$q_{\eta\overline{\eta}\eta}$','$q_{\eta\overline{\eta}\overline{\eta}}$',...
          'phii','phih',...
          'wi','wh','ci','ch',...
          'sigma','siges','sigek',...
          'muri','murh',...
          'muk','muq','r',...
          'muni','munh','signis','signhs','signik','signhk',...
          '$\sigma^{\xi}$','$\sigma^{\xi}$','$sigxik$','$sigxhk$','$\iota^i$','$\iota^h$','$\mu^z$'};

latexv = {'$v^i$','$v_\eta^i$','$v_z^i$','$v^h$','$v_\eta^h$','$v_z^h$'};

%% Initial guess for variables:
elseif(section == "Guess") & strcmp(par.guess,'off')
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
    vi = NaN(par.n1,par.n2);
    vh = NaN(par.n1,par.n2);
    
    for i1=1:par.n1
        for i2=1:par.n2
            vi(i1,i2) = par.rhoi;
            vh(i1,i2) = par.rhoh;
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
    
elseif (section == "Model")
%% Model specific set up:
%state variables in this model are eta (e) and sigma (z):
sigma = z; 

wi = psii/e;
wh = (1-psii)/(1-e);

%foc wrt to consumption:
ci    = vi^((1-zetai)/(1-gammai));
ch    = vh^((1-zetah)/(1-gammah));

%investment maximization:
iotai  = (q-1)/kappa_p;
iotah  = (q-1)/kappa_p;

% per setup of the model:
phii  = log(1+kappa_p*iotai)/kappa_p - deltai;
phih  = log(1+kappa_p*iotah)/kappa_p - deltah;

muz = kappa_z*(sigma - sigbar);

% from aggregation:
muk   = psii*phii + (1-psii)*phih;
    
% Volatilities 
% derived from model set up
signis = wi*sigqs;
signhs = wh*sigqs;

signik = wi*(sigqk+sigma);
signhk = wh*(sigqk+sigma);

siges  = signis -  sigqs       ;
sigek  = signik - (sigqk+sigma);

% value function is a function of the two state variables sigma (z) and eta
% (e)
sigxik = vei/vi*sigek*e;
sigxhk = veh/vh*sigek*e;

sigxis = vei/vi*siges*e + vzi/vi*sigz*z;
sigxhs = veh/vh*siges*e + vzh/vh*sigz*z;

muq   =  qe/q*mue*e + qz/q*muz*z ...
         + 1/2*qee/q*( (siges*e)^2 + (sigek*e)^2 )...
         + 1/2*qzz/q*(sigz*z)^2 ...
         + qez*siges*e*sigz*z;       

muri  = (ai-iotai)/q + phii + muq + sigma*sigqk;
murh  = (ah-iotah)/q + phih + muq + sigma*sigqk;

r = muri ...
    - gammai*wi*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
    + sigqs*sigxis + (sigqk+sigma)*sigxik;

muni  = r + wi*(muri-r) - ci; 
munh  = r + wh*(murh-r) - ch;


% System to Solve with Newton-Raphson
eqmue = ( muni - muk - muq - sigma*sigqk + 1/2*(sigqk+sigma)^2 + 1/2*sigqs^2 - wi*sigqs^2 - wi*(sigqk+sigma)^2 ) - mue; %#ok<*NASGU>

eqq   = (ci*e + ch*(1-e))*q - psii*(ai-iotai) - (1-psii)*(ah-iotah); 

eqpsii = muri - murh... 
        + gammah*wh*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
        - gammai*wi*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
        + sigqs*sigxis + (sigqk+sigma)*sigxik ...
        - sigqs*sigxhs - (sigqk+sigma)*sigxhk;
       
eqsigqs = (sigz*z*qz + siges*e*qe) - sigqs*q;
eqsigqk = sigek*e*qe - sigqk*q;

end


