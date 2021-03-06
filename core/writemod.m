
par.tmpFolder = ['tmp_',par.name];
if not(isfolder(par.tmpFolder)) 
    mkdir(par.tmpFolder)
end

path_tmp = regexp(path,';','split');
if any(find(contains(path_tmp,[pwd,'\tmp_'])))
    itmp = find(contains(path_tmp,[pwd,'\tmp_']));
    for ip = 1:length(itmp)
        rmpath(path_tmp{itmp(ip)});
    end
end

if ispc
    addpath([pwd,'\',par.tmpFolder])
else
    addpath([pwd,'/',par.tmpFolder])
end

if strcmp(par.search,'off')
    modfile = regexp(fileread(['mod_',par.name,'.m']),'\n','split');
elseif strcmp(par.search,'on')
    modfile = regexp(fileread([par.tmpFolder,'/mod_',par.name,'_search.m']),'\n','split');
end

parfile = regexp(fileread(['par_',par.name,'.m']),'\n','split');

linepar = find(contains(modfile,'%% PARAMETERS'));
linevar = find(contains(modfile,'%% VARIABLES'));
linegue = find(contains(modfile,'%% GUESS'));
linemod = find(contains(modfile,'%% MODEL'));
linecon = find(contains(modfile,'%% CONSTRAINTS'));
linesol = find(contains(parfile,'%% SOLVER PARAMETERS'));

cd(par.tmpFolder)
fid = fopen('mod_parameters.m','w');
for i=linepar+1:linevar-1
    fprintf(fid,'%s',modfile{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('mod_variables.m','w');
for i=linevar+1:linegue-1
    fprintf(fid,'%s',modfile{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('mod_guess.m','w');
for i=linegue+1:linemod-1
    fprintf(fid,'%s',modfile{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('mod_model.m','w');
for i=linemod+1:linecon-1
    fprintf(fid,'%s',modfile{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('mod_constraints.m','w');
for i=linecon+1:length(modfile)
    fprintf(fid,'%s',modfile{i});
    fprintf(fid,'\n');
end
fclose(fid);

fid = fopen('par_solver.m','w');
for i=linesol+1:length(parfile)
    fprintf(fid,'%s',parfile{i});
    fprintf(fid,'\n');
end
fclose(fid);

cd ..
