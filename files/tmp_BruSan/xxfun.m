SS  = reshape(S,par.ns,par.N);
DD  = reshape(D,par.nd,par.N);
VV  = reshape(V0,par.nv,par.N);
XX  = reshape(X0,par.nx,par.N);
LL0 = reshape(L0,par.nl,par.N);
MM1 = reshape(M1,par.nl,par.N);
MM2 = reshape(M2,par.nl,par.N);
PP1 = reshape(P1,par.nl,par.N);
PP2 = reshape(P2,par.nl,par.N);

XX_ = Xfun11(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par);
X_ = squeeze(reshape(XX_,[par.nx_ par.dim]));

Xtmp = squeeze(reshape(Xfun00(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,1,1) = Xtmp(:,1,1);

Xtmp = squeeze(reshape(Xfun01(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,1,2:end-1) = Xtmp(:,1,2:end-1);

Xtmp = squeeze(reshape(Xfun02(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,1,end) = Xtmp(:,1,end);

Xtmp = squeeze(reshape(Xfun10(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,2:end-1,1) = Xtmp(:,2:end-1,1);

Xtmp = squeeze(reshape(Xfun12(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,2:end-1,end) = Xtmp(:,2:end-1,end);

Xtmp = squeeze(reshape(Xfun20(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,end,1) = Xtmp(:,end,1);

Xtmp = squeeze(reshape(Xfun21(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,end,2:end-1) = Xtmp(:,end,2:end-1);

Xtmp = squeeze(reshape(Xfun22(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
X_(:,end,end) = Xtmp(:,end,end);

if vtmp==0
    X_ = squeeze(reshape(Xfun33(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));
end

XX_  = reshape(X_,par.nx_,par.N);
