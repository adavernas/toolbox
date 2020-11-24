function J = HJBmap2D(DD,J,r,u,mu,sig,par)
 
dm = DD(1:par.ndim,:);
dp = DD(1+par.ndim:2*par.ndim,:);
da = DD(1+2*par.ndim:3*par.ndim,:);

eta   = cell(1,par.N);
Seta  = NaN(1,par.N);
xi    = cell(1,par.N);
testD = cell(1,par.N);

sig = squeeze(reshape(sig,[par.ns*(par.ns+1)/2 par.dim]));
sig(:,1:par.maxP,:) = 0;
sig(:,end:end-par.maxP,:) = 0;
sig(:,:,1:par.maxP) = 0;
sig(:,:,end:end-par.maxP) = 0;
sig = reshape(sig,par.ns*(par.ns+1)/2,par.N);

Sig = cell(1,par.N);
for ik=1:par.N
    Sig{ik} = NaN(par.ndim,par.ndim);
    for id=1:par.ndim
    for jd=1:par.ndim
        if id==jd
            Sig{ik}(id,jd) = sig(id,ik)/(2*da(id,ik)*da(jd,ik));
        else
            Sig{ik}(id,jd) = sig(3,ik)/(2*da(id,ik)*da(jd,ik));
        end
    end
    end
end

mdim = NaN(1,par.ndim);
for id=1:par.ndim
    mdim(id) = prod([1 par.dim(1:id-1)]);
end

if par.ns>1
    for ik=1:par.N
        [eta{ik},xi{ik},testD{ik}] = decompfun(Sig{ik},par);
        
        for ix=1:length(eta{ik})
            xdim = sum(mdim*xi{ik}(:,ix));
            if or(ik-xdim<1,ik-xdim>par.N)
                eta{ik}(:) = 0;
            end
            if or(ik+xdim<1,ik+xdim>par.N)
                eta{ik}(:) = 0;
            end
        end
        
        Seta(1,ik) = sum(eta{ik});
    end
else
    for ik=1:par.N
        eta{ik} = Sig{ik};
        xi{ik} = 1;
        Seta(1,ik) = Sig{ik};
    end
end
    
Dk =  min(0,mu)./dm;
Uk = -max(0,mu)./dp;

Dk = squeeze(reshape(Dk,[par.ns par.dim]));
Dk(1,1,:) = 0;
if par.ns==2
    Dk(2,:,1) = 0;
end
Dk = reshape(Dk,par.ns,par.N);

Uk = squeeze(reshape(Uk,[par.ns par.dim]));
Uk(1,end,:) = 0;
if par.ns==2
    Uk(2,:,end) = 0;
end
Uk = reshape(Uk,par.ns,par.N);

Mk  = r + 1/par.dt - sum(Dk,1) - sum(Uk,1) + 2*Seta;
Mk_ = r + 1/par.dt - sum(Dk,1) - sum(Uk,1);

indd = ones([par.ns par.dim]);
indd(:,  [1:par.nextr end-par.nextr+1:end],:) = 0;
indd(:,:,[1:par.nextr end-par.nextr+1:end])   = 0;
ind = reshape(indd,par.ns,par.N);
ind = ind(1,:);
% indd = squeeze(indd(1,:,:));

adj1 = zeros(par.dim);
adj1(1:par.nextr        ,:) = repmat( (par.nextr:-1:1)',1,par.dim(2));
adj1(end-par.nextr+1:end,:) = repmat(-(1:par.nextr)'   ,1,par.dim(2));

adj2 = zeros(par.dim);
adj2(:,1:par.nextr)         = repmat( (par.nextr:-1:1),par.dim(1),1);
adj2(:,end-par.nextr+1:end) = repmat(-(1:par.nextr)   ,par.dim(1),1);

A = zeros(par.N,par.N);
for ik=1:par.N
    if or(par.ns==1,ind(ik)==1)
        A(ik,ik) = A(ik,ik)+Mk(ik);
        for id=1:par.ndim
            if ik-mdim(id)>=1
                A(ik,ik-mdim(id)) = A(ik,ik-mdim(id))+Dk(id,ik);
            end
            if ik+mdim(id)<=par.N
                A(ik,ik+mdim(id)) = A(ik,ik+mdim(id))+Uk(id,ik);
            end
        end
        for ix=1:length(eta{ik})
            xdim = sum(mdim*xi{ik}(:,ix));
            if and(ik-xdim>=1,ik-xdim<=par.N)
                A(ik,ik-xdim) = A(ik,ik-xdim)-eta{ik}(ix);
            end
            if and(ik+xdim>=1,ik+xdim<=par.N)
                A(ik,ik+xdim) = A(ik,ik+xdim)-eta{ik}(ix);
            end
        end
        
    elseif ind(ik)==0
    
        inddmap = zeros(1,par.N);
        inddmap(ik) = 1;
        [d1,d2] = find(reshape(inddmap,[par.dim]));
        indmap = zeros(par.dim);
        indmap(d1+adj1(d1,d2),d2+adj2(d1,d2)) = 1;
        ik_ = find(reshape(indmap,1,par.N));
        
        A(ik,ik)  = A(ik,ik)  + Mk_(ik);
        A(ik,ik_) = A(ik,ik_) + 2*Seta(ik_);
        for id=1:par.ndim
            if ik-mdim(id)>=1
                A(ik,ik-mdim(id)) = A(ik,ik-mdim(id))+Dk(id,ik);
            end
            if ik+mdim(id)<=par.N
                A(ik,ik+mdim(id)) = A(ik,ik+mdim(id))+Uk(id,ik);
            end
        end
        for ix=1:length(eta{ik_})
            xdim = sum(mdim*xi{ik_}(:,ix));
            if and(ik_-xdim>=1,ik_-xdim<=par.N)
                A(ik,ik_-xdim) = A(ik,ik_-xdim)-eta{ik_}(ix);
            end
            if and(ik_+xdim>=1,ik_+xdim<=par.N)
                A(ik,ik_+xdim) = A(ik,ik_+xdim)-eta{ik_}(ix);
            end
        end
    end
end

A = sparse(A);
F = A\(u' + 1/par.dt.*J');

% JJt = squeeze(reshape(F',[1 par.dim]));
% JJ  = squeeze(reshape(J,[1 par.dim]));
% 
% rr = squeeze(reshape(r,[1 par.dim]));
% mmu = squeeze(reshape(mu(1,:),[1 par.dim]));
% 
% ddm = squeeze(reshape(dm(1,:),[1 par.dim]));
% ddp = squeeze(reshape(dp(1,:),[1 par.dim]));
% 
% JJp = [JJt(:,2:end) JJt(:,end)];
% JJm = [JJt(:,1) JJt(:,1:end-1)];
% 
% rr.*JJt - ( u + max(mmu,0)./ddp.*(JJp-JJt) + min(mmu,0)./ddm.*(JJt-JJm) + (JJ-JJt)/par.dt )
% 
% [squeeze(reshape(diag(A)+0,[1 par.dim]));rr + 1./par.dt + max(mmu,0)./ddp - min(mmu,0)./ddm]
% 
% A*F - u' - 1/par.dt.*J'

J = F';




