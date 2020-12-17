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

% volatilities derived from model set up
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

% system to Solve with Newton-Raphson
eqmue = ( muni - muk - muq - sigma*sigqk + 1/2*(sigqk+sigma)^2 + 1/2*sigqs^2 - wi*sigqs^2 - wi*(sigqk+sigma)^2 ) - mue; %#ok<*NASGU>

eqq   = (ci*e + ch*(1-e))*q - psii*(ai-iotai) - (1-psii)*(ah-iotah);

eqpsii = muri - murh...
    + gammah*wh*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
    - gammai*wi*( (sigqs)^2 + (sigqk+sigma)^2 ) ...
    + sigqs*sigxis + (sigqk+sigma)*sigxik ...
    - sigqs*sigxhs - (sigqk+sigma)*sigxhk;

eqsigqs = (sigz*z*qz + siges*e*qe) - sigqs*q;
eqsigqk = sigek*e*qe - sigqk*q;

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
     + gammai/2*( wi*(sigqs)^2 + wi*(sigqk+sigma)^2 ) );

rVh   = - (1-gammah)*( 1/(1-1/zetah)*(ch-rhoh) ...
     + r - ch ...
     + gammah/2.*( wh*(sigqs)^2 + wh*(sigqk+sigma)^2 ) );

