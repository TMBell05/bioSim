% bioSimV6

% creates bat-like flocking behavior exiting cave.


% clear all
% close all
% clc

% savename = 'science';
% 
% if ~isempty(dir([savename '*.mat']))
%     files = dir([savename '*.mat']);
%     fprintf('Previous files found.\nBuilding on: %s\n',files(end).name)
%     load(files(end).name,'popNum','caveMouthSz','avgSpeed','cvSpacing','minSpacing','collSpacing',...
%         'dt','timeSteps','viewAng','dispHt','targPt','batLoc','batVel','Pring','Plead','Plost','Pst','roi','maxHts')
%     tmp1           = batLoc(:,1,end);
%     tmp2           = batLoc(:,2,end);
%     tmp3           = batLoc(:,3,end);
%     tmp4           = batVel(:,1,end);
%     tmp5           = batVel(:,2,end);
%     tmp6           = batVel(:,3,end);
%     clear batLoc batVel
%     
%     % preallocate space for arrays
%     batLoc      = nan(popNum,3,timeSteps);
%     batVel      = nan(popNum,3,timeSteps);
%     
%     % initialize bat locations relative to cave mouth (x,y,z) [m]
%     batLoc(:,1,1) = tmp1;
%     batLoc(:,2,1) = tmp2;
%     batLoc(:,3,1) = tmp3;
%     
%     % initialize bat velocities (u,v,w) [m/s]
%     batVel(:,1,1) = tmp4;
%     batVel(:,2,1) = tmp5;
%     batVel(:,3,1) = tmp6;
%     
%     % add in max heights
%     batLoc(:,3,end) = maxHts;
%     
%     clear tmp*
%     files = str2double(files(end).name(length(savename)+1));
%     savename = [savename num2str(files+1)];
%     clear files
% else
%     
%     % simulation settings
%     savename    = [savename '1'];
%     popNum      = 50000;   % size of cave population                    [#]
%     caveMouthSz = [5 3];   % width and height of cave mouth             [m]
%     cvSpacing   = 1;       % minimum flight spacing while exiting cave  [m]
%     minSpacing  = 1;       % minimum flight spacing                     [m]
%     cvSpeed     = 20;       % flight speed while exiting cave            [m/s]
%     avgSpeed    = 30;       % average flight speed                       [m/s]
%     collSpacing = 0.005;   % spacing of in-air collision                [m]
%     dt          = 5;       % update time                                [s]
%     timeSteps   = 100;    % number of simulation time steps            [#]
%     viewAng     = 45;      % cone of view angle                         [degrees]
%     roi         = 100;      % radius of influence                        [m]
%     dispHt      = 2000;    % dispersion height                          [m]
%     targPt      = [2000 -2000 2500]; % target location (x,y,z)           [m]
% 
%     
%     % weighting vector [cohesion alignment seperation original target]
%     Pst         = [1 1 nan .5 .7];
% 
%     % weighting vector for leader class
%     Plead       = [-.5 1 -.1 .9 .7];
%     
%     % weighting vector for lost class
%     Plost       = [1 nan -1 .4 .9];
%     
%     % weighting vector for divergent behavior
%     Pring       = [-1 1 .5 .5 nan];
%     
%     
%     % preallocate space for arrays
%     batLoc      = single(nan(popNum,3,timeSteps));
%     batVel      = single(nan(popNum,3,timeSteps));
%    
%     % maximum number of bats allowed in plane of mouth of cave at once [#]
%     xnum = floor((caveMouthSz(1)-cvSpacing)/cvSpacing*2);
%     ynum = floor((caveMouthSz(2)-cvSpacing)/cvSpacing*2);
%     maxPlaneNum = xnum*ynum;
%     
%     % initialize bat locations relative to cave mouth (x,y,z) [m]
%     runWayLen     = popNum/maxPlaneNum*cvSpacing;
%     
%     tmp           = repmat(linspace(-caveMouthSz(1)/2,caveMouthSz(1)/2,xnum),[1 ceil(popNum/xnum)])';
%     batLoc(:,1,1) = tmp(1:popNum);
%     tmp           = repmat(kron(linspace(-caveMouthSz(2)/2,caveMouthSz(2)/2,ynum),ones(1,xnum))',[ceil(popNum/xnum/ynum) 1]);
%     batLoc(:,2,1) = tmp(1:popNum);
%     tmp           = kron(linspace(0,-runWayLen,ceil(popNum/xnum/ynum)),ones(1,xnum*ynum))';
%     batLoc(:,3,1) = tmp(1:popNum);
%     
%     batLoc(:,1,1) = batLoc(:,1,1)-2*batLoc(:,3,1);
%     batLoc(:,2,1) = batLoc(:,2,1)-2*batLoc(:,3,1);
%     
%     
%     % add slight random component
%     if ~isempty(dir('randSeed.mat'))
%         fprintf('loading randSeed.mat\n')
%         load('randSeed.mat')
%     else
%         fprintf('creating randSeed.mat\n')
%         rdSeed = .5*rand(100000,3);
%         save('randSeed.mat','rdSeed')
%     end    
%     batLoc(:,1,1) = batLoc(:,1,1)+rdSeed(1:size(batLoc,1),1);
%     batLoc(:,2,1) = batLoc(:,2,1)+rdSeed(1:size(batLoc,1),2);
%     batLoc(:,3,1) = batLoc(:,3,1)+rdSeed(1:size(batLoc,1),3);
%     
%     % initialize bat velocities (u,v,w) [m/s]
%     batVel(:,:,1) = -cvSpeed/sqrt(3)*ones(popNum,3);
%     batVel(:,3,1) = -batVel(:,3,1)/2;
%     
%     % initialize maximum height vector
%     maxHts = batLoc(:,3,1);
% end



% for t = 90 : timeSteps
%     disp(t)
%     
%     % update location based on last velocity (x,y,z) from (u,v,w)
%     batLoc(:,1,t) = batLoc(:,1,t-1)+ dt*batVel(:,1,t-1);
%     batLoc(:,2,t) = batLoc(:,2,t-1)+ dt*batVel(:,2,t-1);
%     batLoc(:,3,t) = batLoc(:,3,t-1)+ dt*batVel(:,3,t-1);
%         
%     for batCtr = 1:popNum
%         if max(batLoc(batCtr,3,:))<0 % cave exit
%             
%  
%             % u
%             batVel(batCtr:end,1,t) = batVel(batCtr:end,1,t-1);
%             % v
%             batVel(batCtr:end,2,t) = batVel(batCtr:end,2,t-1);
%             % w
%             batVel(batCtr:end,3,t) = batVel(batCtr:end ,3,t-1);
%             break
%      
%             
%         elseif max(batLoc(batCtr,3,:))<dispHt % emergence column
%             
%             minSpacing = batLoc(batCtr,3,t-1)/1000;
%             
%             % set standard weights
%             P = Pst;
%             
%             % find bats location relative to current bat
%             relX = batLoc(:,1,t)-batLoc(batCtr,1,t);
%             relY = batLoc(:,2,t)-batLoc(batCtr,2,t);
%             relZ = batLoc(:,3,t)-batLoc(batCtr,3,t);
%             currBatRel = sqrt(relX.^2+relY.^2+relZ.^2);
%                         
%             % find vector away from nearest neighbor
%             temp    = sort(currBatRel,'ascend');
%             nearIdx = find(currBatRel==temp(2),1,'first');
%             near(1) = -relX(nearIdx);
%             near(2) = -relY(nearIdx);
%             near(3) = -relZ(nearIdx);
%                         
%             % find normalized seperation velocity vector
%             Vsp(1) = near(1)/norm(near);
%             Vsp(2) = near(2)/norm(near);
%             Vsp(3) = near(3)/norm(near);
%             
%             % check for collision
%             if currBatRel(nearIdx) <= collSpacing
%                 % in-air collision; rebound in opposite direction
%                 batVel(batCtr,1,t) = avgSpeed*Vsp(1);
%                 batVel(batCtr,2,t) = avgSpeed*Vsp(2);
%                 batVel(batCtr,3,t) = avgSpeed*Vsp(3);
%                 continue
%             elseif currBatRel(nearIdx) <= minSpacing && max(batLoc(batCtr,3,:))>0
%                 % getting too close; add seperation
%                 P(3) = 10*(minSpacing-currBatRel(nearIdx))/minSpacing;
%             end
%             
%             % find angles relative to current bat
%             theta = real(acosd((batVel(batCtr,1,t-1)*relX+batVel(batCtr,2,t-1)*relY+batVel(batCtr,3,t-1)*relZ)./...
%                 norm(batVel(batCtr,:,t-1)).*norm([relX relY relZ])));
% 
%             
%             % check who is around the bat (within veiw angle and roi)
%             if sum(currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2))==0
%                
%                 if sum(currBatRel<=roi & currBatRel~=0)==0  % lost. return to group at same layer.
%                    
%                     P      = Plost;
%                     % find center of mass location of bats within radius of influence
%                     com(1) = nanmean(batLoc((batLoc(:,3,t)>batLoc(batCtr,3,t)-100)&(batLoc(:,3,t)<batLoc(batCtr,3,t)+100),1,t)-batLoc(batCtr,1,t));
%                     com(2) = nanmean(batLoc((batLoc(:,3,t)>batLoc(batCtr,3,t)-100)&(batLoc(:,3,t)<batLoc(batCtr,3,t)+100),2,t)-batLoc(batCtr,2,t));
%                     com(3) = nanmean(batLoc((batLoc(:,3,t)>batLoc(batCtr,3,t)-100)&(batLoc(:,3,t)<batLoc(batCtr,3,t)+100),3,t)-batLoc(batCtr,3,t));
%                     Vnear  = [nan nan nan];
%                     
%                 else    % leader. 
%                     
%                     P      = Plead;
%                     % find center of mass location of bats within radius of influence
%                     com(1) = nanmean(batLoc((currBatRel<=roi & currBatRel~=0),1,t)-batLoc(batCtr,1,t));
%                     com(2) = nanmean(batLoc((currBatRel<=roi & currBatRel~=0),2,t)-batLoc(batCtr,2,t));
%                     com(3) = nanmean(batLoc((currBatRel<=roi & currBatRel~=0),3,t)-batLoc(batCtr,3,t));
% %                     com    = [nan nan nan];
%                     
%                     % find average velocity velocity vector of neighbors
%                     Vnear(1) = mean(batVel((currBatRel<=roi & currBatRel~=0),1,t-1));
%                     Vnear(2) = mean(batVel((currBatRel<=roi & currBatRel~=0),2,t-1));
%                     Vnear(3) = mean(batVel((currBatRel<=roi & currBatRel~=0),3,t-1));
%                     
%                 end
%                 
%             else
%                 
%                 % find center of mass location of bats within radius of influence
%                 com(1) = nanmean(batLoc((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),1,t)-batLoc(batCtr,1,t));
%                 com(2) = nanmean(batLoc((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),2,t)-batLoc(batCtr,2,t));
%                 com(3) = nanmean(batLoc((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),3,t)-batLoc(batCtr,3,t));
%                 
%                 % find average velocity velocity vector of neighbors
%                 Vnear(1) = mean(batVel((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),1,t-1));
%                 Vnear(2) = mean(batVel((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),2,t-1));
%                 Vnear(3) = mean(batVel((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),3,t-1));
%                 
%             end
%             
%             % find normalized cohesion velocity vector
%             Vco(1) = com(1)/norm(com);
%             Vco(2) = com(2)/norm(com);
%             Vco(3) = com(3)/norm(com);
%             
%             % find normalized target velocity vector
%             targ(1) = targPt(1) - batLoc(batCtr,1,t);
%             targ(2) = targPt(2) - batLoc(batCtr,2,t);
%             targ(3) = targPt(3) - batLoc(batCtr,3,t);
%             
%             Vtarg(1) = targ(1)/norm(targ);
%             Vtarg(2) = targ(2)/norm(targ);
%             Vtarg(3) = targ(3)/norm(targ);
%             
%             
%             if batLoc(batCtr,3,t-1)< 75
%                 batLoc(batCtr,3,t) = abs(batLoc(batCtr,3,t));
%             end
%             
% 
%         else  % ring spreading
%             
%             P = Pring; 
%               
%             
%             % find bats location relative to current bat
%             relX = batLoc(:,1,t)-batLoc(batCtr,1,t);
%             relY = batLoc(:,2,t)-batLoc(batCtr,2,t);
%             relZ = batLoc(:,3,t)-batLoc(batCtr,3,t);
%             currBatRel = sqrt(relX.^2+relY.^2+relZ.^2);
%             
%             % find angles relative to current bat
%             theta = real(acosd((batVel(batCtr,1,t-1)*relX+batVel(batCtr,2,t-1)*relY+batVel(batCtr,3,t-1)*relZ)./...
%                 norm(batVel(batCtr,:,t-1)).*norm([relX relY relZ])));
%              
%             % find vector away from nearest neighbor
%             temp    = sort(currBatRel,'ascend');
%             nearIdx = find(currBatRel==temp(2),1,'first');
%             near(1) = -relX(nearIdx);
%             near(2) = -relY(nearIdx);
%             near(3) = -relZ(nearIdx);
%             
%             % find normalized seperation velocity vector
%             Vsp(1) = near(1)/norm(near);
%             Vsp(2) = near(2)/norm(near);
%             Vsp(3) = near(3)/norm(near);
% 
%              % check for collision
%             if currBatRel(nearIdx) <= collSpacing
%                 % in-air collision; rebound in opposite direction
%                 batVel(batCtr,1,t) = avgSpeed*Vsp(1);
%                 batVel(batCtr,2,t) = avgSpeed*Vsp(2);
%                 batVel(batCtr,3,t) = avgSpeed*Vsp(3);
%                 continue
%             elseif currBatRel(nearIdx) <= minSpacing && max(batLoc(batCtr,3,:))>0
%                 % getting too close; add seperation
%                 P(3) = 10*(minSpacing-currBatRel(nearIdx))/minSpacing;
%             end
%                 
%              % keep in layer
%              if batLoc(batCtr,3,t) > (dispHt+500*rand(1))
%                  com(1) = nanmean(batLoc(batLoc(:,3,t)>100&batLoc(:,3,t)<(dispHt+100),1,t)-batLoc(batCtr,1,t));
%                  com(2) = nanmean(batLoc(batLoc(:,3,t)>100&batLoc(:,3,t)<(dispHt+100),2,t)-batLoc(batCtr,2,t));
%                  com(3) = batLoc(batCtr,3,t-1)+10;
% 
%              elseif batLoc(batCtr,3,t) < (50*rand(1))
%                  com(1) = nanmean(batLoc(batLoc(:,3,t)>(dispHt-10)&batLoc(:,3,t)<(dispHt+10),1,t)-batLoc(batCtr,1,t));
%                  com(2) = nanmean(batLoc(batLoc(:,3,t)>(dispHt-10)&batLoc(:,3,t)<(dispHt+10),2,t)-batLoc(batCtr,2,t));
%                  com(3) = -norm(com(1:2))/(1+10*rand(1));
% 
%              else
%                  % within layer
%                  com(1) = nanmean(batLoc(batLoc(:,3,t)>(batLoc(batCtr,3,t)-10)&batLoc(:,3,t)<(batLoc(batCtr,3,t)+10),1,t)-batLoc(batCtr,1,t));
%                  com(2) = nanmean(batLoc(batLoc(:,3,t)>(batLoc(batCtr,3,t)-10)&batLoc(:,3,t)<(batLoc(batCtr,3,t)+10),2,t)-batLoc(batCtr,2,t));
%                  com(3) = -rand(1)*norm(com(1:2))*batVel(batCtr,3,t-1)/avgSpeed;
%                  p(4)   = 1;
%              end
%              
%             % find normalized cohesion velocity vector
%             Vco(1) = com(1)/norm(com);
%             Vco(2) = com(2)/norm(com);
%             Vco(3) = com(3)/norm(com);
%             
%             if isnan(Vco)
%                 P(1) = nan;
%             end
%             
%             % find average velocity velocity vector of neighbors
%             Vnear(1) = mean(batVel((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),1,t-1));
%             Vnear(2) = mean(batVel((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),2,t-1));
%             Vnear(3) = mean(batVel((currBatRel<=roi & currBatRel~=0 & (theta<=viewAng/2 | theta>=360-viewAng/2)),3,t-1));
%             
%             Vtarg = [nan nan nan];
%             
%         end
% 
%         
% 
%             
%             % u
%             batVel(batCtr,1,t) = ...
%                 nansum([avgSpeed*P(1)*Vco(1) P(2)*Vnear(1) avgSpeed*P(3)*Vsp(1)...
%                 P(4)*batVel(batCtr,1,t-1) avgSpeed*P(5)*Vtarg(1)])/...
%                 nansum(abs(P));
%             % v
%             batVel(batCtr,2,t) = ...
%                 nansum([avgSpeed*P(1)*Vco(2) P(2)*Vnear(2) avgSpeed*P(3)*Vsp(2)...
%                 P(4)*batVel(batCtr,2,t-1) avgSpeed*P(5)*Vtarg(2)])/...
%                 nansum(abs(P));
%             % w
%             batVel(batCtr,3,t) = ...
%                 nansum([avgSpeed*P(1)*Vco(3) P(2)*Vnear(3) avgSpeed*P(3)*Vsp(3)...
%                 P(4)*batVel(batCtr,3,t-1) avgSpeed*P(5)*Vtarg(3)])/...
%                 nansum(abs(P));
%           
%     end
%     maxHts = max(max(batLoc(:,3,:),[],3),maxHts);
%     
% end




% disp(['saving: ' savename '.mat'])
% save([savename '.mat'],'popNum','caveMouthSz','avgSpeed','cvSpacing','minSpacing','collSpacing',...
%         'dt','timeSteps','viewAng','dispHt','targPt','batLoc','batVel','Pring','Plead','Plost','Pst','roi','maxHts')

% break






close all
scrsz = get(0,'ScreenSize');
cubeSize = 5000; % [m]
ht       = 2000; % [m]

figure('Position',[1 scrsz(4)/1.1 scrsz(3)/1.1 scrsz(4)/1.1])
for ct = 1:1000
    
    subplot(1,3,1) 
    plot(batLoc(:,1,ct),batLoc(:,2,ct),'k.')
    xlim([-cubeSize/2 cubeSize/2])
    ylim([-cubeSize/2 cubeSize/2])
    title('x/y reference frame')
    grid on
    axis square
    
    subplot(1,3,2)
    plot(batLoc(:,1,ct),batLoc(:,3,ct),'k.')
    xlim([-cubeSize/2 cubeSize/2])
    ylim([-10 ht])
    title('x/z  reference frame')
    grid on
    
    subplot(1,3,3)
    plot3(batLoc(:,1,ct),batLoc(:,2,ct),batLoc(:,3,ct),'k.')
    xlim([-cubeSize/2 cubeSize/2])
    ylim([-cubeSize/2 cubeSize/2])
    zlim([0 ht])
    title(num2str(ct))
    grid on
    shg
  pause(.5)
end






close all
scrsz = get(0,'ScreenSize');
cubeSize = 6; % [km]
ht       = 1750; % [m]

%figure('Position',[1 scrsz(4)/4 scrsz(3)/4 scrsz(4)/2])

nn = '1200';
%for ct = 1:10:1000
    ct = 200;
    figure
    plot(batLoc(:,1,ct)/1000,batLoc(:,2,ct)/1000,'k.')
    xlim([-cubeSize/2 cubeSize/2])
    ylim([-cubeSize/2 cubeSize/2])
    title('Horizontal Projection of Agent Locations','fontsize',24)
    xlabel('longitudinal distance from cave mouth [km] ','fontsize',18)
    ylabel('latitudinal distance from cave mouth [km]','fontsize',18)
    set(gca,'XTick',-10:10)
    set(gca,'YTick',-10:10)
    axis square
    grid on
    saveas(gcf,['/Users/phillipstepanian/research/documents/h' nn '.eps'])
    
    
    figure
    plot(batLoc(:,1,ct)/1000,batLoc(:,3,ct),'k.')
    xlim([-cubeSize/2 cubeSize/2])
    ylim([0 ht])
    title('West-East Vertical Projection of Agent Locations','fontsize',24)
    xlabel('longitudinal distance from cave mouth [km] ','fontsize',18)
    ylabel('height above ground level [m]  ','fontsize',18)
    set(gca,'XTick',-10:10)
    set(gca,'YTick',0:250:2750)
    axis square
    grid on
    saveas(gcf,['/Users/phillipstepanian/research/documents/v' nn '.eps'])
   
   
%     
%     subplot(1,3,3)
%     plot3(batLoc(:,1,ct),batLoc(:,2,ct),batLoc(:,3,ct),'k.')
%     xlim([-cubeSize/2 cubeSize/2])
%     ylim([-cubeSize/2 cubeSize/2])
%     zlim([0 ht])
%     title(num2str(ct))
%     grid on
    shg
     %pause(.5)
 
%end








ct = 50;

figure(1)
colormap(boonlib('zmapn'))
plot(batLoc(:,1,ct)/1000+bOrX/1000,batLoc(:,2,ct)/1000+bOrY/1000,'k.','MarkerSize',5)
axis([-5 70 -10 65])
%axis('square')
title({'Biological Simulation ';'Horizontal Projection of Agent Locations  '},'fontsize',16)
xlabel('longitudinal distance from radar [km]  ','fontsize',14)
ylabel('latitudinal distance from radar [km]  ','fontsize',14)
hold on
plot(0,0,'o','color',[0 0 0],'MarkerSize',10)



% 
% 
% 
% 
% 
% 
% for t = 2 : timeSteps
%     disp(t)
%     
%     % update location based on last velocity (x,y,z) from (u,v,w)
%     batLoc(:,1,t) = batLoc(:,1,t-1)+ dt*batVel(:,1,t-1);
%     batLoc(:,2,t) = batLoc(:,2,t-1)+ dt*batVel(:,2,t-1);
%     batLoc(:,3,t) = batLoc(:,3,t-1)+ dt*batVel(:,3,t-1);
%     
%     
%     batVel(:,:,t) = batVel(:,:,t-1);
%     
%     nm = sqrt(batVel(:,1,t).^2+batVel(:,2,t).^2+batVel(:,3,t).^2);
%     batVel(nm<3,1,t)=2*batVel(nm<3,1,t);
%     batVel(nm<3,2,t)=2*batVel(nm<3,2,t);
%     
%     tp = 1000+1400*rand(100000,1);
%     bt = 50*rand(100000,1);
%     
%     batVel(batLoc(:,3,t)<bt & batVel(:,3,t)<0,3,t)=-0.9*batVel(batLoc(:,3,t)<bt & batVel(:,3,t)<0,3,t-1);
%     batVel(batLoc(:,3,t)<bt & batVel(:,3,t)<0,1,t)=1.1*batVel(batLoc(:,3,t)<bt & batVel(:,3,t)<0,1,t-1);
%     batVel(batLoc(:,3,t)<bt & batVel(:,3,t)<0,2,t)=1.1*batVel(batLoc(:,3,t)<bt & batVel(:,3,t)<0,2,t-1);
%     
%     
%     batVel(batLoc(:,3,t)>tp & batVel(:,3,t-1)>0,3,t)=-batVel(batLoc(:,3,t)>tp & batVel(:,3,t)>0,3,t-1); 
%     
% end
