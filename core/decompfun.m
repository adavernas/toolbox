function [eta,xi,testD,B] = decompfun(A,par)

n = size(A,1);
if any(any(tril(A,-1)~=triu(A,1)'))
    disp('Not a triangular diffusion matrix!')
end

neta = n + 2*(n^2-n)/2;
eta  = NaN(1,neta);
xi   = zeros(n,neta);

[~,I] = sort(diag(A),1,'descend');
Ah = A;
for i=1:2
    for j=1:2
        Ah(i,j) = A(I(i),I(j));
    end
end
ncorr = 0;
if Ah(2,1)<0
    Ah = abs(Ah);
    ncorr = 1;
end
    
e = 0;
for i=1:n
    e = e+1;
    A_ = Ah;
    A_(i,i) = 0;
    eta(e)  = Ah(i,i) - sum(abs(A_(i,:)));
    xi(i,e) = 1;
end

for l=1:n
for i=1:l-1
    e = e+1;
    eta(e)  = Ah(i,l)*(Ah(i,l)>=0);
    xi(i,e) = 1;
    xi(l,e) = 1;
end
end

for l=1:n
for i=1:l-1
    e = e+1;
    eta(e)  = -Ah(i,l)*(Ah(i,l)<0);
    xi(i,e) =  1;
    xi(l,e) = -1;
end
end
    
if  and(~isdiagdom(Ah),n<=2)
   
    q  = 0;
    p  = 1;
    qp = 1;
    pp = 1;

    v_ah = vecfun(Ah);
%     n_ah = norm(v_ah);
%     v_ah = v_ah/n_ah;
    
    iterD = 0;
    stop  = 0;
    while stop==0
        iterD = iterD+1;
        
        xxi  = [p; q];
        xxip = [pp;qp];
        X  = xxi *xxi';
        Xp = xxip*xxip';
         
        v_x  = vecfun(X);
%         n_x  = norm(v_x);
%         v_x  = v_x/n_x;
        
        v_xp = vecfun(Xp);
%         n_xp = norm(v_xp);
%         v_xp = v_xp/n_xp;

        v_n   = cross(v_x,v_xp); % normal vector
        v_ap  = Pfun(v_n,zeros(3,1),v_ah); % projection on plane
        
        if or(p+pp>par.maxP,norm(v_ah-v_ap)<=par.tolD*norm(v_ah))
                        
            eta = [v_x v_xp]\v_ap;
            xi = [xxi xxip];
            stop = 1;
            
%              eta(1)*v_x + eta(2)*v_xp - vecfun(Ah)
%              eta(1)*X + eta(2)*Xp - Ah
%             
%              eta(1)*v_x + eta(2)*v_xp - v_ap*n_ah
%              eta(1)*X/n_x + eta(2)*Xp/n_xp - matfun(v_ap*n_ah)
%              
             
%             eta(1)*v_x + eta(2)*v_xp - v_ap
            
%         elseif if(iterD>par.maxD)
% %             
%             testQ = 1;
%             iterQ = 0;
%             while and(iterQ<par.maxQ,testQ>par.tolD)
%                 iterQ = iterQ+1;
%                 
%                 xi  = permn(-par.maxP:par.maxP,n)';
%                 nxi = size(xi,2);
%                 xi  = xi(:,randi(nxi,1,par.dimQ));
%                 
%                 options = optimoptions(@fmincon,'Algorithm','interior-point','TolFun',par.tolD,'StepTolerance',par.tolD);
%                 options = optimoptions(options,'Display','off');
%                 
%                 eta0 = zeros(1,par.dimQ);
%                 eta = fmincon(@(eta) decfun(eta,A,xi),eta0,[],[],[],[],zeros(1,par.dimQ),[],[],options);
%                 stop = 1;
%                 testQ = decfun(eta,A,xi);
%             end
            
        else
            qpp = q+qp;
            ppp = p+pp;
            xxipp = [ppp;qpp];
            Xpp = xxipp*xxipp';
            v_xpp = vecfun(Xpp);
            
            eta_ = [v_x v_xp v_xpp]\v_ah;
            
            v_np = cross(v_x,v_xpp); % normal vector
            s = sign(dot(v_np,[1/2;0;1/2])); % direction of normal vector
                
            if all(eta_>=0)
  
                eta = eta_;
                xi = [xxi xxip xxipp];
                stop = 1;
                
            elseif s*dot(v_np,v_ah)<=0 % if not in same half-space as point (1,0,1)
                qp = qpp;
                pp = ppp;
            else
                q = qpp;
                p = ppp;
            end
            
%             plottest
%             drawnow
        end
    end
end   
xi = xi(I,:);
if ncorr
    xi(2,1) = xi(2,1)*-1; 
    xi(2,2) = xi(2,2)*-1; 
    if length(eta)==3
    	xi(2,3) = xi(2,3)*-1; 
    end
end

[testD,~,~,B] = decfun(eta,A,xi);

if testD>1
    disp(['stencil error: ',num2str(testD)])
end
% disp(eta')
  
%   'ok'