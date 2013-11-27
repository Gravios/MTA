function [walk_state,walk_feature] = walk(Session,method,varargin)

switch method

  case 'vel'
    [walkThresh,minimum_interval,display] = DefaultArgs(varargin,{2.5,60,0});
    % Walk - speed threshold
    % Feature - spine_lower 
    fet = xyzvel(Session.xyz,gmi(Session.Model.markers,'spine_lower'),[1 2])*Session.xyzSampleRate/10;
    walk_feature = medfilt1(fet,20);
    walk = zeros(size(walk_feature));
    walk(walk_feature>walkThresh) = 1;
    walkPeriods = ThreshCross(walk,0.5,minimum_interval);

  case 'com0415'
    [walkThresh,minimum_interval,display] = DefaultArgs(varargin,{0.0316,60,0});
    if isempty(Session.ang),
        Session = Session.load_ang;
    end
    win = 7;
    xyz = reshape(Filter0(gausswin(win)./sum(gausswin(win)),Session.xyz),size(Session.xyz,1),size(Session.xyz,2),size(Session.xyz,3));
    v = sqrt(sum(diff(xyz).^2,3));
    win = 241;
    comb = Filter0(gausswin(win)./sum(gausswin(win)),Session.com(Session.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'})));
    comh = Filter0(gausswin(win)./sum(gausswin(win)),Session.com(Session.Model.rb({'head_back','head_left','head_front','head_right'})));
    comb = sqrt(sum(diff(comb(:,[1,2])).^2,2));
    comh = sqrt(sum(diff(comh(:,[1,2])).^2,2));
    comv = [comb,comh];
    for i =1:2,
        comv(:,i) = ButFilter(comv(:,i),11,2./(Session.xyzSampleRate/2),'low');
        comv(:,i) = clip(comv(:,i),0.003,30);
    end
    sfet = [circ_dist(Session.ang(:,2,3,1),Session.ang(:,1,2,1)),...
            circ_dist(Session.ang(:,3,4,1),Session.ang(:,2,3,1)),...
            circ_dist(Session.ang(:,4,5,1),Session.ang(:,3,4,1)),...
            circ_dist(Session.ang(:,5,7,1),Session.ang(:,4,5,1))];
    num_fet = size(sfet,2);
    for i = 1:num_fet
        sfet(isnan(sfet(:,i)),i) = circ_mean(sfet(~isnan(sfet(:,i)),i));
    end
    sfet = Filter0(gausswin(7)./sum(gausswin(7)),sfet);
    bsfet = ButFilter(sum(sfet,2),11,6./(Session.xyzSampleRate/2),'low');
    tsfet = Filter0(gausswin(141)./sum(gausswin(141)),bsfet);
    dbsfet = diff(bsfet-tsfet); 
    nfet = -Filter0(gausswin(141)./sum(gausswin(141)),abs(dbsfet.*Filter0(gausswin(141)./sum(gausswin(141)),comv(:,1)).*12))./(Session.xyz(1:end-1,7,3)./Session.ang(1:end-1,3,4,2)).*1000;
    walk_feature = cat(1,nfet(1),nfet);
    walk = zeros(size(walk_feature));
    walk(walk_feature>walkThresh) =1;
    walkPeriods = ThreshCross(walk,0.5,minimum_interval);    

  case 'com'
    [walkThresh,minimum_interval,display] = DefaultArgs(varargin,{2,40,0});
    % Walk - speed threshold
    % Feature - center of mass of the body markers
    b4 = Session.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
    fet = Session.com(b4);

    fetv = sqrt(sum(diff(fet(:,:,[1 2]),1).^2,3));
    walk_feature = medfilt1(fetv,32).*Session.xyzSampleRate./10;
    walk = zeros(size(walk_feature));
    walk(walk_feature>walkThresh) =1;
    walkPeriods = ThreshCross(walk,0.5,minimum_interval);

  case 'head'
    [walkThresh,minimum_interval,display] = DefaultArgs(varargin,{0.6,64,0});
    % Walk - speed threshold
    % Feature - center of mass of the head markers
    b4 = Session.Model.rb({'head_back','head_left','head_front','head_right'});
    fet = Session.com(b4);

    fetv = sqrt(sum(diff(fet(:,:,[1 2]),1).^2,3));
    walk_feature = medfilt1(fetv,32);
    walk = zeros(size(walk_feature));
    walk(walk_feature>walkThresh) =1;
    walkPeriods = ThreshCross(walk,0.5,minimum_interval);

  case 'fet'
    b4 = Session.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
    fet = Session.com(b4);
    fetv = sqrt(sum(diff(fet(:,:,[1 2]),1).^2,3));
    walkPeriods = medfilt1(fetv,32).*Session.xyzSampleRate./10;
    
    

  otherwise
    error(['Method: ' method ' does not exist'])
end


walk_state = walkPeriods;