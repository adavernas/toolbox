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

