function [rV,uV,muV,sigV] = HJB(SS,VV,XX,XX_,par,j)

e = SS(par.S.ie,:);
z = SS(par.S.iz,:);

mue  = XX(par.X.imue,:);
sigqs = XX(par.X.isigqs,:);
sigqk = XX(par.X.isigqk,:);

muz  = XX_(par.X_.imuz,:);

siges = XX_(par.X_.isiges,:);
sigek = XX_(par.X_.isigek,:);

r     = XX_(par.X_.ir,:);
sigma = XX_(par.X_.isigma,:);

sigz = par.sigz;


if strcmp(j,'i')
    rho    = par.rhoi;
    zeta   = par.zetai;
    gamma  = par.gammai;
    
    c      = XX_(par.X_.ici,:);
    w      = XX_(par.X_.iwi,:);
    murj   = XX_(par.X_.imuri,:);

    v      = VV(par.V.ivi,:);
    ve     = VV(par.V.ivei,:);
    vz     = VV(par.V.ivzi,:);
    
elseif strcmp(j,'h')
    rho    = par.rhoh;
    zeta   = par.zetah;
    gamma  = par.gammah;
    
    c      = XX_(par.X_.ich,:);
    w      = XX_(par.X_.iwh,:);
    murj   = XX_(par.X_.imurh,:);
    
    v      = VV(par.V.ivh,:);
    ve     = VV(par.V.iveh,:);  
    vz     = VV(par.V.ivzh,:);  
end



uV   = zeros(1,par.N);
muV  = [mue.*e;muz.*z];
sigV = [(siges.*e).^2 + (sigek.*e).^2;
        (sigz.*z).^2;
         siges.*e.*sigz.*z];


   
rV   = - (1-gamma)*( 1./(1-1/zeta).*(c-rho) ...
     + r - c ...
     + gamma./2.*( w.*(sigqs).^2 + w.*(sigqk+sigma).^2 ) );




