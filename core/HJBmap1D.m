function Jp = HJBmap1D(DD,J,r,u,mu,sig,par)
 
dm = reform(DD(1:par.ndim,:),par);
dp = reform(DD(1+par.ndim:2*par.ndim,:),par);

J   = reform(J,par);
r   = reform(r,par);
u   = reform(u,par); 
mu  = reform(mu,par); 
sig = reform(sig,par);

U = cell(par.ndim^2,1);
D = cell(par.ndim^2,1);
M = cell(par.ndim^2,1);
A = cell(par.ndim^2,1);
F = cell(par.ndim^2,1);
for is=1:par.ndim
for ib=1:par.ndim
    
in = (is-1)*par.ndim+ib;

% sig{in}([1 end]) = 0;
U{in} = - max(mu{in},0)./dp{in} - sig{in}.^2.*dm{in}./(dm{in}.*dp{in}.*(dm{in}+dp{in}));
U{in}(end) = 0;

D{in} =   min(mu{in},0)./dm{in} - sig{in}.^2.*dp{in}./(dm{in}.*dp{in}.*(dm{in}+dp{in}));
D{in}(1) = 0;

M{in}  = r{in} + 1/par.dt - U{in} - D{in};
% M{in}([1 end]) = 0;

A{in} = spdiags(M{in}',0,par.dim(is),par.dim(is)) ...
      + spdiags([0 U{in}(1:end-1)]',1,par.dim(is),par.dim(is)) ...
      + spdiags([D{in}(2:end) 0]',-1,par.dim(is),par.dim(is));

% A{in}(end,end)   = 1;
% A{in}(end,end-1) = -1;
% 
% A{in}(1,1) = 1;
% A{in}(1,2) = -1;
% 
% J{in}([1 end]) = 0;
% u{in}([1 end]) = 0;

% R*J = u + mu*Jx + sig^2*Jxx/2 + Jt
% A*J = u + J/dt

F{in} = NaN(size(A{in}));
if all(~isnan(A{in}))
    F{in} = A{in}\(u{in}' + 1/par.dt.*J{in}');
end

end
end

Jp = NaN(par.dim);
Jp(~isnan(F{1}),1)   = F{1}(~isnan(F{1}));
Jp(~isnan(F{2}),end) = F{2}(~isnan(F{2}));
Jp(1,~isnan(F{3}))   = F{3}(~isnan(F{3}));
Jp(end,~isnan(F{4})) = F{4}(~isnan(F{4}));

% vec1 = linspace(0.01,0.99,par.n1); % ok lineraire
% vec1 = 3*vec1.^2 - 2*vec1.^3;
% vec2 = linspace(0.01,.99,par.n2);
% vec2 = 3*vec2.^2 - 2*vec2.^3;
% 
% [x1, x2]  = ndgrid(vec1, vec2);
% 
% ncol = 2;
% nlin = 2;
% 
% figure(5); clf(5); h = 0;
% 
% h = h+1;
% subplot(ncol,nlin,h); hold on
% Mp = NaN(par.dim);
% Mp(:,1)   = M{1};
% Mp(:,end) = M{2};
% Mp(1,:)   = M{3};
% Mp(end,:) = M{4};
% m = mesh(x1,x2,Mp);
% set(m,'facecolor','none')
% view([-32.4,42])
% xlabel('\eta')
% ylabel('z');  xlim([vec1(1) vec1(end)]); ylim([vec2(1) vec2(end)]);
% title('$q_\eta$','Interpreter','LateX');
% 
% h = h+1;
% subplot(ncol,nlin,h); hold on
% mup = NaN(par.dim);
% mup(:,1)   = mu{1};
% mup(:,end) = mu{2};
% mup(1,:)   = mu{3};
% mup(end,:) = mu{4};
% m = mesh(x1,x2,mup);
% set(m,'facecolor','none')
% view([-32.4,42])
% xlabel('\eta')
% ylabel('z');  xlim([vec1(1) vec1(end)]); ylim([vec2(1) vec2(end)]);
% title('$q_\eta$','Interpreter','LateX');

