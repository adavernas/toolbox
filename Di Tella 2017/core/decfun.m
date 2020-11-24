function [F,dF,H,B] = decfun(eta,A,xi)
% Compute the error, the gradient and Hessian of the stencil decomposition

n = size(xi,1);
neta = length(eta);
eb = NaN(n,n,neta);
b  = NaN(n,n,neta);
for e=1:neta
    eb(:,:,e) = eta(e)*xi(:,e)*xi(:,e)';
    b(:,:,e)  = xi(:,e)*xi(:,e)';
end
B = sum(eb,3);
F = max(max( (A - B) ))./sqrt(A(1,1)^2+A(2,2)^2);

dF = NaN(neta,1);
for e=1:neta
    dF(e) = -sum(sum( (A - B).*b(:,:,e) ));
end

H = NaN(neta,neta);
for e1=1:neta
for e2=1:neta
    H(e1,e2) = sum(sum( b(:,:,e1).*b(:,:,e2) ));
end
end
