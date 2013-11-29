

Trial = MTATrial('jg05-20120317');
fwin = gausswin(11)./sum(gausswin(11));
sst = 'w';
Trial.ang.load(Trial.sync);
Trial.xyz.filter(fwin);


vel =[];
vel = [vel,Filter0(gausswin(61)./sum(gausswin(61)),log10(Trial.vel({'spine_lower'}))).^2];
vel = [vel,log10(Trial.vel({'spine_lower','spine_upper','head_front'},[1,2]))];
vel = [vel,log10(Trial.vel({'spine_lower','spine_upper','head_front'},3))];
vel = [vel,sqrt(sum(diff(Trial.com(Trial.xyz.model.rb({'spine_lower','pelvis_root','spine_middle'}))).^2,3))];
%vel = [vel,Filter0(fwin,diff(circ_dist(Trial.ang(1:end,1,2,1),Trial.ang(1:end,1,3,1)))+diff(circ_dist(Trial.ang(1:end,2,3,1),Trial.ang(1:end,2,4,1)))+diff(circ_dist(Trial.ang(1:end,3,4,1),Trial.ang(1:end,3,5,1))))];
vel = [vel,log10(abs(Filter0(fwin,circ_dist(Trial.ang(1:end-1,1,2,1),Trial.ang(1:end-1,1,3,1)))).*...
          abs(Filter0(fwin,circ_dist(Trial.ang(1:end-1,2,3,1),Trial.ang(1:end-1,2,4,1)))).*...
          abs(Filter0(fwin,circ_dist(Trial.ang(1:end-1,3,4,1),Trial.ang(1:end-1,3,5,1)))))];

% aba = abs(Filter0(fwin,circ_dist(Trial.ang(1:end-1,1,2,1),Trial.ang(1:end-1,1,3,1)))).*...
%       abs(Filter0(fwin,circ_dist(Trial.ang(1:end-1,2,3,1),Trial.ang(1:end-1,2,4,1)))).*...
%       abs(Filter0(fwin,circ_dist(Trial.ang(1:end-1,3,4,1),Trial.ang(1:end-1,3,5,1))));
% 
% sfet = [circ_dist(Trial.ang(1:end-1,2,3,1),Trial.ang(1:end-1,1,2,1)),...
%         circ_dist(Trial.ang(1:end-1,3,4,1),Trial.ang(1:end-1,2,3,1)),...
%         circ_dist(Trial.ang(1:end-1,4,5,1),Trial.ang(1:end-1,3,4,1)),...
%         circ_dist(Trial.ang(1:end-1,5,7,1),Trial.ang(1:end-1,4,5,1))];
       
if isempty(y),
dav = Filter0(fwin,diff(circ_dist(Trial.ang(1:end,1,2,1),Trial.ang(1:end,1,3,1))))+...
      Filter0(fwin,diff(circ_dist(Trial.ang(1:end,2,3,1),Trial.ang(1:end,2,4,1))))+...
      Filter0(fwin,diff(circ_dist(Trial.ang(1:end,3,4,1),Trial.ang(1:end,3,5,1))));
wdav = WhitenSignal(dav);
[y,t,f] = mtchglong(wdav,2^8,Trial.xyz.sampleRate,2^7,2^7-1,[],[],[],[1,20]);
end
%vel = [vel,[log10(sum(y(:,f>8&f<16),2)./sum(y(:,f<7),2));zeros(Trial.xyz.size(1)-numel(t)-1,1)]];
  
vel = [vel,Trial.ang(1:end-1,1,2,2)];
vel = [vel,Trial.ang(1:end-1,3,4,2)];
vel = [vel,vel.^2];

vel = MTADxyz([],[],Filter0(fwin,vel),Trial.xyz.sampleRate);


vel_walk = vel(Trial.stc{'w'}(1:100,:),:);
vel_walk(sum(isnan(vel_walk)|isinf(vel_walk),2)>0,:) = [];
vel_rear = vel(Trial.stc{'r'}(1:50,:),:);
vel_rear(sum(isnan(vel_rear)|isinf(vel_rear),2)>0,:) = [];



cov_walk = cov(vel_walk);
cov_rear = cov(vel_rear);

vel.data(vel.data==0) = nan;
vel_not_nan = ~isnan(vel.data(:,1));

mean_vel_walk = repmat(nanmean(vel_walk),vel.size(1),1);
d_walk = zeros(vel.size(1),1);
d_walk(vel_not_nan) = -.5*log(det(cov_walk))-.5*dot(((vel.data(vel_not_nan,:)-mean_vel_walk(vel_not_nan,:))/cov_walk)',(vel.data(vel_not_nan,:)-mean_vel_walk(vel_not_nan,:))');

mean_vel_rear = repmat(nanmean(vel_rear),vel.size(1),1);
d_rear = zeros(vel.size(1),1);
d_rear(vel_not_nan) = -.5*log(det(cov_walk))-.5*dot(((vel.data(vel_not_nan,:)-mean_vel_rear(vel_not_nan,:))/cov_rear)',(vel.data(vel_not_nan,:)-mean_vel_rear(vel_not_nan,:))');

dwin= gausswin(21)./sum(gausswin(21));
figure
%plot(Filter0(dwin,d_walk),'.'),hold on,plot(Filter0(dwin,d_rear),'.r'),
plot(Filter0(dwin,d_walk)),hold on,plot(Filter0(dwin,d_rear),'r'),
Lines(Trial.stc{'w',Trial.xyz.sampleRate}(:),[],'k');
Lines(Trial.stc{'r',Trial.xyz.sampleRate}(:),[],'r');
Lines([],0,'k')
