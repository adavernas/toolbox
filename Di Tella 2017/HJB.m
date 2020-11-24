function [rV,uV,muV,sigV] = HJB(SS,VV,XX,XX_,par,j)

e = SS(par.S.ie,:);
z = SS(par.S.iz,:);

mue  = XX(par.X.imue,:);
sigq = XX_(par.X_.isigq,:);

muz  = XX_(par.X_.imuz,:);


sigx = XX(par.X.isigx,:);
sigz_ = XX_(par.X_.isigz_,:);

r     = XX_(par.X_.ir,:);
pii = XX_(par.X_.ipii,:);

%sigz = par.sigz;
psii = par.psii ; 
gamma  = par.gammai;


if strcmp(j,'i')
    rho    = par.rhoi;
    zeta   = par.zetai;
    
    c      = XX_(par.X_.ici,:);
    sigxij = XX_(par.X_.isigxii,:);
    sigN   = XX_(par.X_.isign,:);
    tet    = XX(par.X.itheta, :);

    v      = VV(par.V.ivi,:);
    ve     = VV(par.V.ivei,:);
    vz     = VV(par.V.ivzi,:);
    
elseif strcmp(j,'h')
    rho    = par.rhoh;
    zeta   = par.zetah;
    
    c      = XX_(par.X_.ich,:);
    sigxij = XX_(par.X_.isigxih,:);

    sigN   = XX(par.X.isigw,:);
    tet    = 0;


    v      = VV(par.V.ivh,:);
    ve     = VV(par.V.iveh,:);  
    vz     = VV(par.V.ivzh,:);  
end



uV   = zeros(1,par.N);
muV  = [mue.*e;muz.*z];
sigV = [(sigx.*e).^2;
        (sigz_.*z).^2;
         sigx.*e.*sigz_.*z];


   
rV   =  - ( 1./(1-psii).*(c - rho) ...
     + r - c  - tet.^2 ...
     + gamma./2.*( (sigN).^2 )) ;




