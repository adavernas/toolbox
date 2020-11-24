function f = clearnan(x,k)

sx = size(x);
f = NaN(sx);

for i=1:sx(1)
for j=1:sx(2)
    if isnan(x(i,j))
        f(i,j) = k;
    elseif ~isnan(x(i,j))
        f(i,j) = x(i,j);
    end
end
end