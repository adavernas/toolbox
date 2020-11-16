function [Jp,DJm,DJp] = finitediff2D(J,D,x1)

D = squeeze(D);
if size(D,2)==1
    D = D';
end

DJp_ = (J(2:end)-J(1:end-1))./D(1:end-1);
DJp = interp1(x1(2:end-1),DJp_(2:end),x1,'spline');

DJm_ = (J(2:end)-J(1:end-1))./D(2:end);
DJm = interp1(x1(2:end-1),DJm_(1:end-1),x1,'spline');

% Jp = DJm.*(mu<0) + DJp.*(mu>0) + 1/2*DJm.*(mu==0) + 1/2*DJp.*(mu==0);
Jp = 1/2*DJm + 1/2*DJp;
