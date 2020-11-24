function isdom = isdiagdom(A)

isdom = true;
for r = 1:size(A,1)
    rowdom = 2*abs(A(r,r)) > sum(abs(A(r,:)));
    isdom = isdom && rowdom;
end