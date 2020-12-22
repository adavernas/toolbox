X0 = NaN([par.nx par.dim]);

X0(par.X.iq,:,:) = 1;
X0(par.X.ipsii,:,:) = 0.95;
X0(par.X.isigqs,:,:) = 0;
X0(par.X.isigqk,:,:) = 0;
X0(par.X.imue,:,:) = 0;

vi = par.rhoi*ones(par.dim);
vh = par.rhoh*ones(par.dim);

