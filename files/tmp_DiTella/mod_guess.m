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

