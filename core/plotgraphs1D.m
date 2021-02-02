aa1 = 1;
bb1 = par.n1;

nlin = floor(sqrt(par.nx));
ncol = ceil((par.nx)/nlin);

figure(1); clf(1);
for h=1:par.nx
    subplot(nlin,ncol,h); hold on
    
    eval(['Zq  = squeeze(X0(par.X.i',vars{h},',:,:));'])
    
    plot(vec1,Zq); grid on;
    
    xlabel(latexs{1},'Interpreter','LateX');

    title(['$',latex{h},'$'],'Interpreter','LateX');
    set(gca,'FontSize',16)
    xlim([vec1(aa1) vec1(bb1)]);
    tix=get(gca,'ztick')';
    set(gca,'zticklabel',num2str(tix,'%.2f'));
end

nlin = floor(sqrt(par.nx_));
ncol = ceil((par.nx_)/nlin);

figure(2); clf(2);
for h_=1:par.nx_
    subplot(nlin,ncol,h_); hold on
    
    eval(['Zq  = squeeze(X_(par.X_.i',vars_{h_},',:,:));'])
    
    plot(vec1,Zq); grid on;
    
    xlabel(latexs{1},'Interpreter','LateX');

    title(['$',latex_{h_},'$'],'Interpreter','LateX');
    set(gca,'FontSize',16)
    xlim([vec1(aa1) vec1(bb1)]);
    tix=get(gca,'ztick')';
    set(gca,'zticklabel',num2str(tix,'%.2f'));
end

nlin = floor(sqrt(par.nv+2));
ncol = ceil((par.nv+2)/nlin);

% value  = {'vi','vh','vi_','vh_','q_','q0'};

ylvd = [0 0 0 0 0.5 0.5];
ylvu = [1000 3 2 2 1.5 1.5];

figure(3); clf(3);
for v=1:par.nv
    subplot(nlin,ncol,v); hold on
    
    eval(['Zq  = squeeze(V0(par.V.i',value{v},',:,:));'])
    
    plot(vec1,Zq); grid on;
    xlim([vec1(aa1) vec1(bb1)]);
%     ylim([ylvd(v) ylvu(v)]);
    
    xlabel(latexs{1},'Interpreter','LateX');

    title(['$',latexv{v},'$'],'Interpreter','LateX');
    set(gca,'FontSize',16)
    
    tix=get(gca,'ztick')';
    set(gca,'zticklabel',num2str(tix,'%.2f'));
end

if exist('rrvi')
v = v+1;
subplot(nlin,ncol,v); hold on
Zq  = squeeze(squeeze(reshape(rrvi,[1 par.dim]))) - squeeze(squeeze(reshape(uuvi,[1 par.dim])))./vvi;

plot(vec1,Zq); grid on;

xlabel(latexs{1},'Interpreter','LateX');

title(['$\mu^\theta$'],'Interpreter','LateX');
set(gca,'FontSize',16)
xlim([vec1(aa1) vec1(bb1)]);
tix=get(gca,'ztick')';
set(gca,'zticklabel',num2str(tix,'%.2f'));
end

% v = v+1;
% subplot(nlin,ncol,v); hold on
% nb_ = squeeze(X_(par.X_.inb_,:,:));
% nh_ = squeeze(X_(par.X_.inh_,:,:));
% 
% plot(vec1,nb_,'r',vec1,nh_,'b'); grid on;
% 
% xlabel('$\eta$','Interpreter','LateX');
% 
% % title(['$\mu^\theta$'],'Interpreter','LateX');
% set(gca,'FontSize',16)
% xlim([vec1(aa1) vec1(bb1)]);
% tix=get(gca,'ztick')';
% set(gca,'zticklabel',num2str(tix,'%.2f'));
% 
% eta_ = squeeze(X_(par.X_.ieta_,:,:))
% nb_.*vec1./(nb_.*vec1+nh_.*(1-vec1)) - eta_
% 
% 
% nb_.*vec1+nh_.*(1-vec1) - q_.*squeeze(X_(par.X_.iEomega1,:,:))

% if strcmp(par.savegraph_,'on')
%     set(gcf,'PaperPosition', [0 0 12 12]);
%     set(gcf,'PaperSize', [12 12]);    
%     filename = [par.name,'_theta_',num2str(par.theta,'%0.0f')];
%     
%     cdtmp = cd;
%     cd ..
%     cd ..
%     graphpath = [pwd,'/graphs/'];
%     cd(cdtmp)
%     
%     saveas(gcf,[graphpath,filename],'pdf')
% end


