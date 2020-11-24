function [Jp,DJm,DJp] = finitediff2D(J,D,x1,x2,d)

D = squeeze(D);
if size(D,1)==1
    D = D';
end

if d==1
    DJp_ = (J(2:end,:)-J(1:end-1,:))./D(1:end-1,:);
    DJp = [DJp_;DJp_(end,:)];
    
%     DJp = interp2(x2(2:end-1,:),x1(2:end-1,:),DJp_(2:end,:),x2,x1,'linear');
    
    DJm_ = (J(2:end,:)-J(1:end-1,:))./D(2:end,:);
    DJm = [DJm_(1,:);DJm_];
    
%     DJm = interp2(x2(2:end-1,:),x1(2:end-1,:),DJm_(1:end-1,:),x2,x1,'linear');

elseif d==2
    DJp_ = (J(:,2:end)-J(:,1:end-1))./D(:,1:end-1);
    DJp = [DJp_ DJp_(:,end)];

%     DJp = interp2(x2(:,2:end-1),x1(:,2:end-1),DJp_(:,2:end),x2,x1,'linear');
        
    DJm_ = (J(:,2:end)-J(:,1:end-1))./D(:,2:end);
    DJm = [DJm_(:,1) DJm_];

%     DJm = interp2(x2(:,2:end-1),x1(:,2:end-1),DJm_(:,1:end-1),x2,x1,'linear');
end

% Jp = DJm.*(mu<0) + DJp.*(mu>0) + 1/2*DJm.*(mu==0) + 1/2*DJp.*(mu==0);
Jp = 1/2*DJm + 1/2*DJp;
