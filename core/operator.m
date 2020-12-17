function O = operator(eta,xi,k,V)
nx = length(eta);

O = 0;
for ix=1:nx
    O = O + eta(ix)*(V(k+xi(:,ix)) + V(k-xi(:,ix)));
end