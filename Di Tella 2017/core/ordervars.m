function fun_ = ordervars(fun,names,short,param,state,dvec,value,vars,varl0,varm1,varm2,varp1,varp2,fvars1,type) %#ok<INUSL,STOUT>

if strcmp(type,'{i}')
    if not(isempty(fvars1))
        vec = char(fvars1(1));
        for i=2:length(fvars1)
            vec = [vec,',',char(fvars1(i))]; %#ok<AGROW>
        end
        
        for j=1:length(names)
            for i=1:length(eval(names{j}))
                vec = strrep(vec,[',',eval([names{1,j},'{i}']),','],[',',short{1,j},'{',num2str(i),'}(ind),']);
            end
            
            for i=1:length(eval(names{j}))
                istr = strfind(vec,',');
                if strcmp(vec(1:istr(1)-1),eval([names{1,j},'{i}']))
                    vec = [short{1,j},'{',num2str(i),'}(ind),',vec(istr(1)+1:end)];
                end
                if strcmp(vec(istr(end)+1:end),eval([names{1,j},'{i}']))
                    vec = [vec(1:istr(end)-1),',',short{1,j},'{',num2str(i),'}(ind)'];
                end
            end
        end
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par,ind)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_(',vec,');';]);
        fclose(fid);
        cd('..')
        
    else
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par,ind)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_;';]);
        fclose(fid);
        cd('..')
    end
    
elseif strcmp(type,'(i,:)')
    if not(isempty(fvars1))
        vec = char(fvars1(1));
        for i=2:length(fvars1)
            vec = [vec,',',char(fvars1(i))]; %#ok<AGROW>
        end
        
        for j=1:length(names)
            for i=1:length(eval(names{j}))
                vec = strrep(vec,[',',eval([names{1,j},'{i}']),','],[',',short{1,j},'(',num2str(i),',:),']);
            end
            
            for i=1:length(eval(names{j}))
                istr = strfind(vec,',');
                if strcmp(vec(1:istr(1)-1),eval([names{1,j},'{i}']))
                    vec = [short{1,j},'(',num2str(i),',:),',vec(istr(1)+1:end)];
                end
                if strcmp(vec(istr(end)+1:end),eval([names{1,j},'{i}']))
                    vec = [vec(1:istr(end)-1),',',short{1,j},'(',num2str(i),',:)'];
                end
            end
        end
        
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_(',vec,');';]);
        fclose(fid);
        cd('..')
        
    else
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_;';]);
        fclose(fid);
        cd('..')
    end
    
elseif strcmp(type,'(i,:,:)')
    if not(isempty(fvars1))
        vec = char(fvars1(1));
        for i=2:length(fvars1)
            vec = [vec,',',char(fvars1(i))]; %#ok<AGROW>
        end
        
        for j=1:length(names)
            for i=1:length(eval(names{j}))
                vec = strrep(vec,[',',eval([names{1,j},'{i}']),','],[',',short{1,j},'(',num2str(i),',:,:),']);
            end
            
            for i=1:length(eval(names{j}))
                istr = strfind(vec,',');
                if strcmp(vec(1:istr(1)-1),eval([names{1,j},'{i}']))
                    vec = [short{1,j},'(',num2str(i),',:,:),',vec(istr(1)+1:end)];
                end
                if strcmp(vec(istr(end)+1:end),eval([names{1,j},'{i}']))
                    vec = [vec(1:istr(end)-1),',',short{1,j},'(',num2str(i),',:,:)'];
                end
            end
        end
        
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_(',vec,');';]);
        fclose(fid);
        cd('..')
        
    else
        
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_;';]);
        fclose(fid);
        cd('..')
    end
    
elseif strcmp(type,'(i)')
    if not(isempty(fvars1))
        vec = char(fvars1(1));
        for i=2:length(fvars1)
            vec = [vec,',',char(fvars1(i))]; %#ok<AGROW>
        end
        
        for j=1:length(names)
            for i=1:length(eval(names{j}))
                vec = strrep(vec,[',',eval([names{1,j},'{i}']),','],[',',short{1,j},'(',num2str(i),'),']);
            end
            
            for i=1:length(eval(names{j}))
                if ~isempty(strfind(['(',vec,')'],['(',eval([names{1,j},'{i}']),')']))
                    vec = strrep(vec,eval([names{1,j},'{i}']),[short{1,j},'(',num2str(i),')']);
                end
            end
            
            for i=1:length(eval(names{j}))
                istr = strfind(vec,',');
                if ~isempty(istr)
                    if strcmp(vec(1:istr(1)-1),eval([names{1,j},'{i}']))
                        vec = [short{1,j},'(',num2str(i),'),',vec(istr(1)+1:end)];
                    end
                    if strcmp(vec(istr(end)+1:end),eval([names{1,j},'{i}']))
                        vec = [vec(1:istr(end)-1),',',short{1,j},'(',num2str(i),')'];
                    end
                end
            end
        end
        
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_(',vec,');';]);
        fclose(fid);
        cd('..')
        
    else
        cd('tmp')
        fid = fopen([fun,'.m'],'w');
        fprintf(fid,'%s',['function F = ',fun,'(P,S,D,V,X,L0,M1,M2,P1,P2,Vone,Vzero,par)']);fprintf(fid,'\n');
        fprintf(fid,'%s',['F = ',fun,'_;';]);
        fclose(fid);
        cd('..')
    end
end


