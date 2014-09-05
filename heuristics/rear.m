function [rear_state] = rear(Session,method,varargin)

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end

ang = Session.ang.copy;
if ang.isempty,
    ang.load(Session);
end

xyz = Session.xyz.copy;
if xyz.isempty,
    xyz.load(Session);
end

if xyz.sampleRate>120,xyz.resample(120);end
if ang.sampleRate>120,ang.resample(120);end

rear_feature = abs(xyz(:,'head_front',3)-xyz(:,'spine_lower',3)).*ang(:,'spine_middle','spine_upper',2);
rear_feature(isnan(rear_feature))=0;
switch method
    
    
    case 'hmm'
        %% Hidden Markov Model
        [minimum_interval,min_inter_rear_duration,hmm_states,display] = DefaultArgs(varargin,{50,10,2,0});
        
        [State, hmm, decode] = gausshmm(rear_feature,hmm_states);
        [~, rearstate ] = max(cat(1,hmm.state.Mu));
        
        rearing_points = zeros(size(rear_feature));
        rearing_points(State==rearstate) = 1;
        rearing_points(1) = 0;
        rearing_points(end) = 0;
        rearPeriods = ThreshCross(rearing_points,0.5,minimum_interval);
        interRearPeriods = [[0,rearPeriods(1,1)];[rearPeriods(1:end-1,2),rearPeriods(2:end,1)];[rearPeriods(end,2),length(rearing_points)]];
        interRearDuration = diff(interRearPeriods,1,2);
        false_rear_off_ind = find(interRearDuration<min_inter_rear_duration);
        if size(false_rear_off_ind,1)>0,
            rpt = rearPeriods;
            rpt(false_rear_off_ind,:) = [rearPeriods(false_rear_off_ind,1) rearPeriods(false_rear_off_ind+1,2)];
            rpt(false_rear_off_ind+1,:) =0;
            rearPeriods = reshape(rpt(find(rpt)),[],2);
        end
        
        
    case 'com'
        %% Can't remember wat COM stands for
        [rearThresh,minimum_interval,display] = DefaultArgs(varargin,{50,64,0});
        rear_feature = MTADxyz('data',rear_feature,'sampleRate',xyz.sampleRate);
        rearPeriods = ThreshCross(rear_feature.filter(gtwin(1,xyz.sampleRate))>rearThresh,...
            0.5,minimum_interval);
        
    case 'com0415'
      %%Non-Operational Code
      [rearThresh,minimum_interval,display] = DefaultArgs(varargin,{0,50,0});
        win = 7;
        xyz = reshape(Filter0(gausswin(win)./sum(gausswin(win)),xyz.data),xyz.size);
        v = sqrt(sum(diff(xyz).^2,3));
        win = 241;
        comb = Filter0(gausswin(win)./sum(gausswin(win)),Session.com(Session.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'})));
        comh = Filter0(gausswin(win)./sum(gausswin(win)),Session.com(Session.model.rb({'head_back','head_left','head_front','head_right'})));
        comb = sqrt(sum(diff(comb(:,[1,2])).^2,2));
        comh = sqrt(sum(diff(comh(:,[1,2])).^2,2));
        comv = [comb,comh];
        for i =1:2,
            comv(:,i) = ButFilter(comv(:,i),11,2./(Session.xyzSampleRate/2),'low');
            comv(:,i) = clip(comv(:,i),0.003,30);
        end
        sfet = [circ_dist(ang(:,2,3,1),ang(:,1,2,1)),...
            circ_dist(ang(:,3,4,1),ang(:,2,3,1)),...
            circ_dist(ang(:,4,5,1),ang(:,3,4,1)),...
            circ_dist(ang(:,5,7,1),ang(:,4,5,1))];
        num_fet = size(sfet,2);
        for i = 1:num_fet
            sfet(isnan(sfet(:,i)),i) = circ_mean(sfet(~isnan(sfet(:,i)),i));
        end
        sfet = Filter0(gausswin(7)./sum(gausswin(7)),sfet);
        bsfet = ButFilter(sum(sfet,2),11,6./(Session.xyzSampleRate/2),'low');
        tsfet = Filter0(gausswin(141)./sum(gausswin(141)),bsfet);
        dbsfet = diff(bsfet-tsfet);
        nfet = -Filter0(gausswin(141)./sum(gausswin(141)),abs(dbsfet.*Filter0(gausswin(141)./sum(gausswin(141)),comv(:,1)).*12))./(Session.xyz(1:end-1,7,3)./ang(1:end-1,3,4,2)).*1000;
        rear_feature = -cat(1,nfet(1),nfet);
        rearing_points = zeros(size(rear_feature));
        rearing_points(rear_feature>rearThresh)=1;
        rearPeriods = ThreshCross(rearing_points,0.5,minimum_interval);
        
    case 'fet'
        rearPeriods = rear_feature;
        
end

rear_state = rearPeriods;


%% Diagnostics
% $$$
% $$$ rint = zeros(size(rearPeriods,1),1);
% $$$ rmax = zeros(size(rearPeriods,1),1);
% $$$ rhmax = zeros(size(rearPeriods,1),1);
% $$$ for i = 1:size(rearPeriods,1),
% $$$ rint(i) = sum(rear_feature(rearPeriods(i,1):rearPeriods(i,2)));
% $$$ rmax(i) = max(rear_feature(rearPeriods(i,1):rearPeriods(i,2)));
% $$$ rhmax(i) = max(Session.xyz(rearPeriods(i,1):rearPeriods(i,2),Session.model.gmi('head_front'),3));
% $$$ end
% $$$ figure
% $$$ hist(log10(rint),100)
% $$$ plot3(log10(rint),rmax,rhmax,'.')
% $$$
% $$$ rear_state = rearPeriods(find(rmax>160),:);
% $$$
% $$$ plot(Session.xyz(:,7,3)),Lines(rear_state(:,1),[],'g');,Lines(rear_state(:,2),[],'r');
% $$$
% $$$ Session.Bhv.States{1}.state = rear_state;
% $$$ Session.Bhv.save(Session,1)
% $$$ Session.save

