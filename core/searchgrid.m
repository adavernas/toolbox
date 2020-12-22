function [] = searchgrid(name,param,grid,varargin)
% This file computes the solution method of d'Avernas and Vandeweyer (2018)
% written by Adrien d'Avernas
% 2019 d'Avernas and Vandeweyer all rights reserved

set(0,'DefaultFigureWindowStyle','docked')
warning('off','MATLAB:interp2:NaNstrip') 

p = inputParser;
addRequired(p,'name',@ischar);
addRequired(p,'param',@ischar);
parse(p,name,param);

par = p.Results;

modfile = regexp(fileread(['mod_',name,'.m']),'\n','split');
line_par_beg = find(contains(modfile,'%% PARAMETERS'));
line_par_end = find(contains(modfile,'%% VARIABLES'));
linepar = find(contains(modfile(line_par_beg:line_par_end),['par.',param,' =']));

par.tmpFolder = ['tmp_',name];
if ~exist(par.tmpFolder, 'dir')
    mkdir(par.tmpFolder)
end

figure(3); clf(3);
for jj=1:length(grid)
    cd(par.tmpFolder)
    fid = fopen(['mod_',par.name,'_search.m'],'w');
    for ii=1:length(modfile)
        if ii==linepar
            fprintf(fid,'%s\n',['par.',param,' = ',num2str(grid(jj)),';']);
        else
            fprintf(fid,'%s\n',modfile{ii});
        end
    end
    fclose(fid);
    cd ..
    
    varargin = [varargin(:)',{'search'},{'on'}];
    
    solvemod(name,varargin{:})
    
    plot_Orchad
end

% searchgrid('Orchad','gammai',[1.2 1.3],'dimensions','1D','guess','on','loop','off','write','off','debug','on','outerplot','on')