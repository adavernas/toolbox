par.ns  = length(state);

if par.ns==2
    names  = {'param','state','dvec','value','vars','varl0','varm1','varm2','varp1','varp2'};
    names_ = {names{:},'vars_','cstv'}; %#ok<CCAT>
    short  = {'P','S','D','V','X','L0','M1','M2','P1','P2'};
    short_ = {short{:},'X_','CC'}; %#ok<CCAT>
elseif par.ns==1
    names  = {'param','state','dvec','value','vars','varl0','varm1','varp1'};
    names_ = {names{:},'vars_','cstv'}; %#ok<CCAT>
    short  = {'P','S','D','V','X','L0','M1','P1'};
    short_ = {short{:},'X_','CC'}; %#ok<CCAT>
end
    
all = [state,value,vars,vars_];
par.np  = length(param);
par.nd  = length(dvec);
par.nv  = length(value);
par.nx  = length(vars);
par.nx_ = length(vars_);
par.nl  = length(last);
par.ncc = length(cstv);
par.nb  = 3^par.ns+1;
par.vc  = 2^par.ncc;
par.nvc = 2^(par.ns + par.ncc);
par.Nk   = prod(par.dim);

varm1 = cell(1,par.nl);
varm2 = cell(1,par.nl);

for il=1:par.nl
	eval(['varl0{il} = ''',last{il},'l0'';']);
end
    
for is=1:par.ns
    for il=1:par.nl
        eval(['varm',num2str(is),'{il} = ''',last{il},'m',state{is},''';']);
    end
end

for is=1:par.ns
    for il=1:par.nl
        eval(['varp',num2str(is),'{il} = ''',last{il},'p',state{is},''';']);
    end
end

for i=1:length(last)
    eval(['par.L.i',char(last{i}),' = ',num2str(i),';']);
end
    
for j=1:length(names_)
    for i=1:length(eval(names_{j}))
        eval(['par.',short_{1,j},'.i',char(eval([names_{1,j},'{i}'])),' = ',num2str(i),';']);
    end
end

if strcmp(par.write,'on')
    fun  = {'dF','F','X','L','CC'};
    fun_ = {'(i)','(i)','(i,:)','(i)','(i)'};
    
    bnd = cell(1,par.nb);
    for ib=1:par.nb-1
        bnd{ib} = dec2base(ib-1,par.base,par.ns);
    end
    
    if par.ns==1
        bnd{end} = '3';
    elseif par.ns==2
        bnd{end} = '33';
    end
    
    if par.ns==1
    vec  = {',names,short,param,state,dvec,value,vars,varl0,varm1,NaN,varp1,NaN,'};
    
    for ibnd = 1:par.nb
        eval(['[dFvars',bnd{ibnd},...
            ',Fvars',bnd{ibnd},...
            ',Xvars',bnd{ibnd},...
            ',Lvars',bnd{ibnd},...
            ',CCvars',bnd{ibnd},...
            '] = writefun(param,state,dvec,value,vars,vars_,last,varl0,varm1,NaN,varp1,NaN,',...
            '''',bnd{ibnd},''',cstv,par);'])
    end
    elseif par.ns==2
    vec  = {',names,short,param,state,dvec,value,vars,varl0,varm1,varm2,varp1,varp2,'};
    
    for ibnd = 1:par.nb
        eval(['[dFvars',bnd{ibnd},...
            ',Fvars',bnd{ibnd},...
            ',Xvars',bnd{ibnd},...
            ',Lvars',bnd{ibnd},...
            ',CCvars',bnd{ibnd},...
            '] = writefun(param,state,dvec,value,vars,vars_,last,varl0,varm1,varm2,varp1,varp2,',...
            '''',bnd{ibnd},''',cstv,par);'])
    end
    end
    
    
    for in = 1:length(fun)
        for ib = 1:par.nb
            eval(['ordervars(''',fun{in},'fun',bnd{ib},'''',vec{:},fun{in},'vars',bnd{ib},',''',fun_{in},''',par);'])
        end
    end
end

% Write Function process1_.m
cd(par.tmpFolder)
if strcmp(par.write,'on')
    fid = fopen('process1_.m','w');
    for ib=1:par.nb
        if ib==1
            fprintf(fid,'%s',['if vtmp==',num2str(base2dec(bnd{ib},par.base)+1)]);fprintf(fid,'\n');
        elseif ib<par.nb
            fprintf(fid,'%s',['elseif vtmp==',num2str(base2dec(bnd{ib},par.base)+1)]);fprintf(fid,'\n');
        elseif ib==par.nb
            fprintf(fid,'%s','elseif vtmp==0');fprintf(fid,'\n');
        end
        
        fprintf(fid,'%s',['    myfun = @(X) Ffun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,real(X),L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);']);fprintf(fid,'\n');        
        fprintf(fid,'%s',['    dFF   = dFfun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);']);fprintf(fid,'\n');
        fprintf(fid,'%s',['    FF    =  Ffun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);']);fprintf(fid,'\n');
        if par.ncc>0
            fprintf(fid,'%s',['    CCtmp = CCfun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);']);fprintf(fid,'\n');
        end
%         if par.nco>0
%             fprintf(fid,'%s',['    COtmp = COfun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);']);fprintf(fid,'\n');
%         end
        if ib==par.nb
            fprintf(fid,'%s','end');fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end
cd('..')

% Write Function process2_.m
if strcmp(par.write,'on')
    cd(par.tmpFolder)
    fid = fopen('process2_.m','w');
    for ib=1:par.nb
        if ib==1
            fprintf(fid,'%s',['if vtmp==',num2str(base2dec(bnd{ib},par.base)+1)]);fprintf(fid,'\n');
        elseif ib<par.nb
            fprintf(fid,'%s',['elseif vtmp==',num2str(base2dec(bnd{ib},par.base)+1)]);fprintf(fid,'\n');
        elseif ib==par.nb
            fprintf(fid,'%s','elseif vtmp==0');fprintf(fid,'\n');
        end
        fprintf(fid,'%s',['    L(:,i1,i2) = Lfun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,XX0,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par);']);
        fprintf(fid,'\n');
        if ib==par.nb
            fprintf(fid,'%s','end');fprintf(fid,'\n');
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    cd('..')
end

% Write Function xxfun.m
if strcmp(par.write,'on')
    cd(par.tmpFolder)
    fid = fopen('xxfun.m','w');
    
    fprintf(fid,'%s',['SS  = reshape(S,par.ns,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['DD  = reshape(D,par.nd,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['VV  = reshape(V0,par.nv,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['XX  = reshape(X0,par.nx,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['LL0 = reshape(L0,par.nl,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['MM1 = reshape(M1,par.nl,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['MM2 = reshape(M2,par.nl,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['PP1 = reshape(P1,par.nl,par.N);']);fprintf(fid,'\n');
    fprintf(fid,'%s',['PP2 = reshape(P2,par.nl,par.N);']);fprintf(fid,'\n');fprintf(fid,'\n');
    if par.ns==1
        fprintf(fid,'%s',['XX_ = Xfun1(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par);']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_ = squeeze(reshape(XX_,[par.nx_ par.dim]));']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun0(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,1) = Xtmp(:,1);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun2(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,end) = Xtmp(:,end);']);fprintf(fid,'\n');fprintf(fid,'\n');

        fprintf(fid,'%s',['if vtmp==0']);fprintf(fid,'\n');
        fprintf(fid,'%s',['    X_ = squeeze(reshape(Xfun3(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['end']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['XX_  = reshape(X_,par.nx_,par.N);']);fprintf(fid,'\n');
        
    elseif par.ns==2
        fprintf(fid,'%s',['XX_ = Xfun11(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par);']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_ = squeeze(reshape(XX_,[par.nx_ par.dim]));']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun00(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,1,1) = Xtmp(:,1,1);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun01(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,1,2:end-1) = Xtmp(:,1,2:end-1);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun02(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,1,end) = Xtmp(:,1,end);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun10(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,2:end-1,1) = Xtmp(:,2:end-1,1);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun12(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,2:end-1,end) = Xtmp(:,2:end-1,end);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun20(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,end,1) = Xtmp(:,end,1);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun21(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,end,2:end-1) = Xtmp(:,end,2:end-1);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['Xtmp = squeeze(reshape(Xfun22(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['X_(:,end,end) = Xtmp(:,end,end);']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['if vtmp==0']);fprintf(fid,'\n');
        fprintf(fid,'%s',['    X_ = squeeze(reshape(Xfun33(P,SS,DD,VV,XX,LL0,MM1,MM2,PP1,PP2,ones(1,par.N),zeros(1,par.N),par),[par.nx_ par.dim]));']);fprintf(fid,'\n');
        fprintf(fid,'%s',['end']);fprintf(fid,'\n');fprintf(fid,'\n');
        
        fprintf(fid,'%s',['XX_  = reshape(X_,par.nx_,par.N);']);fprintf(fid,'\n');
        
    end
    
    fclose(fid);
    cd('..')
end

% Write Function process1fsolve_.m
% options = optimset('display','off','tolF',1e-15);
% if strcmp(par.write,'on')
%     cd(par.tmpFolder)
%     fid = fopen('process1fsolve.m','w');
%     for ib=1:par.nb
%         if ib==1
%             fprintf(fid,'%s',['if vtmp==',num2str(base2dec(bnd{ib},par.base)+1)]);fprintf(fid,'\n');
%         elseif ib<par.nb
%             fprintf(fid,'%s',['elseif vtmp==',num2str(base2dec(bnd{ib},par.base)+1)]);fprintf(fid,'\n');
%         elseif ib==par.nb
%             fprintf(fid,'%s','elseif vtmp==0');fprintf(fid,'\n');
%         end
%         fprintf(fid,'%s',['    X1  =  fmincon(@(XX)sum(Ffun',bnd{ib},'(P,Stmp,Dtmp,V0tmp,XX,L0tmp,M1tmp,M2tmp,P1tmp,P2tmp,1,0,par).^2),XX0,A,B_,[],[],[],[],[],options);']);
%         fprintf(fid,'\n');
%         if ib==par.nb
%             fprintf(fid,'%s','end');fprintf(fid,'\n');
%         end
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
%     cd('..')
% end

% Write Function verfun.m
if strcmp(par.write_con,'on')
    cd(par.tmpFolder)
    fid = fopen('verfun.m','w');
    fprintf(fid,'%s','function [ver,cst1,BC,BC_,icst,icst_] = verfun(SS,XX,CC,cstv,csts,cstn,bndv,bndn,i1,i2,par) %#ok<INUSL>');fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'%s','e  = SS(1);');fprintf(fid,'\n');
    fprintf(fid,'\n');
    if par.ns>1
        fprintf(fid,'%s','z  = SS(2);');fprintf(fid,'\n');
        fprintf(fid,'\n');
        
        fprintf(fid,'%s','if i1==1');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd1 = ''0'';');fprintf(fid,'\n');
        fprintf(fid,'%s','elseif i1==par.n1');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd1 = ''2'';');fprintf(fid,'\n');
        fprintf(fid,'%s','else');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd1 = ''1'';');fprintf(fid,'\n');
        fprintf(fid,'%s','end');fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'%s','if i2==1');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd2 = ''0'';');fprintf(fid,'\n');
        fprintf(fid,'%s','elseif i2==par.n2');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd2 = ''2'';');fprintf(fid,'\n');
        fprintf(fid,'%s','else');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd2 = ''1'';');fprintf(fid,'\n');
        fprintf(fid,'%s','end');fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'%s','bnd = [bnd1,bnd2];');fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'%s','ver = base2dec(bnd,par.base)+1;');
        fprintf(fid,'\n');
    elseif par.ns==1
        fprintf(fid,'%s','if i1==1');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd = ''0'';');fprintf(fid,'\n');
        fprintf(fid,'%s','elseif i1==par.n1');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd = ''2'';');fprintf(fid,'\n');
        fprintf(fid,'%s','else');fprintf(fid,'\n');
        fprintf(fid,'%s','    bnd = ''1'';');fprintf(fid,'\n');
        fprintf(fid,'%s','end');fprintf(fid,'\n');
        fprintf(fid,'\n');
        fprintf(fid,'%s','ver = base2dec(bnd,par.base)+1;');
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'%s','cst1  = NaN(par.ncc,1);');fprintf(fid,'\n');
    fprintf(fid,'\n');
    for icc=1:par.ncc
        if strcmp(csts{icc},'max')
            eq1 = '<';
            eq2 = '>';
        elseif strcmp(csts{icc},'min')
            eq1 = '>';
            eq2 = '<';
        end
        
        fprintf(fid,'%s',['cst1(',num2str(icc),') = not(XX(par.X.i',cstv{icc},')',eq1,cstn{icc},');']);fprintf(fid,'\n');
        fprintf(fid,'\n');
    end
    
    % Setup Constaint Control
    fprintf(fid,'%s','icst  = NaN(par.ncc,1);');
    fprintf(fid,'\n');
    fprintf(fid,'%s','icst_ = NaN(par.ncc,1);');
    fprintf(fid,'\n');
    fprintf(fid,'%s','BC    = NaN(par.ncc,1);');
    fprintf(fid,'\n');
    fprintf(fid,'%s','BC_   = NaN(par.ncc,1);');
    fprintf(fid,'\n');

    for icc=1:par.ncc
    fprintf(fid,'\n');
    fprintf(fid,'%s','icst(',num2str(icc),')   = par.X.i',cstv{icc},';');
    fprintf(fid,'\n');
    fprintf(fid,'%s','icst_(',num2str(icc),')  = par.X.i',cstv_{icc},';');
    fprintf(fid,'\n');
    fprintf(fid,'%s','BC(',num2str(icc),')   = ',cstn{icc},';');
    fprintf(fid,'\n');
    fprintf(fid,'%s','BC_(',num2str(icc),')  = ',cstn_{icc},';');
    fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    
    fclose(fid);

    cd('..')
end

%% State Variables
dp1 = vec1(2:par.n1) - vec1(1:par.n1-1); % diff upward
dp1 = [dp1 dp1(end)];
dm1 = vec1(2:par.n1) - vec1(1:par.n1-1); % diff downward
dm1 = [dm1(1) dm1];

tmp = vec1(2:par.n1) - vec1(1:par.n1-1); % average over -par.maxP to +par.maxP
da1 = NaN(1,par.n1);
for i=1+par.maxP:par.n1-par.maxP
    da1(i) = mean(tmp(i-par.maxP:i+par.maxP-1));
end
da1(1:par.maxP) = da1(1+par.maxP);
da1(par.n1-par.maxP+1:par.n1) = da1(par.n1-par.maxP);

if par.n2>1
    dp2 = vec2(2:par.n2) - vec2(1:par.n2-1);
    dp2 = [dp2 dp2(end)];
    dm2 = vec2(2:par.n2) - vec2(1:par.n2-1);
    dm2 = [dm2(1) dm2];
    
    tmp = vec2(2:par.n2) - vec2(1:par.n2-1);
    da2 = NaN(1,par.n2);
    for i=1+par.maxP:par.n2-par.maxP
        da2(i) = mean(tmp(i-par.maxP:i+par.maxP-1));
    end
    da2(1:par.maxP) = da2(1+par.maxP);
    da2(par.n2-par.maxP+1:par.n2) = da2(par.n2-par.maxP);
    da2(isnan(da2)) = 1/2*dm2(isnan(da2))+1/2*dp2(isnan(da2));
    
    [dm1,dm2] = ndgrid(dm1, dm2);
    [dp1,dp2] = ndgrid(dp1, dp2);
    [da1,da2] = ndgrid(da1, da2);
end

vec1_ = vec1([1+par.nextr:end-par.nextr]);

if par.n2>1
    vec2_ = vec2([1+par.nextr:end-par.nextr]);
    [x1, x2]  = ndgrid(vec1, vec2);
    [x1_,x2_] = ndgrid(vec1_,vec2_);
end

S = NaN([par.ns par.dim]);
if par.n2>1
    S(1,:,:) = x1;
    S(2,:,:) = x2;
else
    S(1,:) = vec1;
end

% CO = NaN([par.nco par.dim]);
C1 = NaN([par.ncc par.dim]);

L  = NaN([par.nl par.dim]);
M1 = NaN([par.nl par.dim]);
M2 = NaN([par.nl par.dim]);
P1 = NaN([par.nl par.dim]);
P2 = NaN([par.nl par.dim]);

X_ = NaN([par.nx_ par.dim]);

D = NaN([par.nd par.dim]);
if par.ns==1
    D(1,:,:) = dm1;
    D(2,:,:) = dp1;
    D(3,:,:) = da1;
    
elseif par.ns==2
    D(1,:,:) = dm1;
    D(2,:,:) = dm2;
    
    D(3,:,:) = dp1;
    D(4,:,:) = dp2;
    
    D(5,:,:) = da1;
    D(6,:,:) = da2;
end

P = NaN(par.np,1);
for ip=1:par.np
    eval(['P(ip) = par.',param{ip},';'])
end
