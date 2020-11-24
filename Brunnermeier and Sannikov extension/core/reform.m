function F = reform(f,par)

ndim = size(f,1);
Ftmp = reshape(f,[ndim par.dim]);

F = cell(4,1);

F{1} = squeeze(Ftmp(1,:,1));
F{2} = squeeze(Ftmp(1,:,end));

if ndim>1
    F{3} = squeeze(Ftmp(2,1,:))';
    F{4} = squeeze(Ftmp(2,end,:))';
else
    F{3} = squeeze(Ftmp(1,1,:))';
    F{4} = squeeze(Ftmp(1,end,:))';
end

