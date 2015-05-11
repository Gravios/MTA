function fet = fet_lgr(Trial,varargin)
%function fet = fet_lgr(Trial,sampleRate)
%


if ~isempty(varargin)
sampleRate = varargin{1};
else 
    sampleRate = Trial.xyz.sampleRate;
end

dsx = Trial.load('xyz');

vl = dsx.vel([1,3,7],[1,2]);
vl.filter('ButFilter',3,4,'low');
vl.data(vl.data<0.0001) = 0.0001;
vl.resample(sampleRate);
vl.data = log10(vl.data);

dsx.resample(sampleRate);
dsx.filter('ButFilter',3,4);
dsa = create(Trial.ang.copy,Trial,dsx);


% turning fet smoothed angular speed in the xy plane
w = round(90/120*sampleRate);
lw = round(70/120*sampleRate);
h = round(w/2);
lh = round(lw/2);
d = circshift(circ_dist(circshift(dsa(:,'pelvis_root','spine_upper',1),-w),dsa(:,'pelvis_root','spine_upper',1)),h);
d = circshift(nanmean(GetSegs(d,1:size(d,1),lw,nan))',lh);
ds = circshift(circ_dist(circshift(dsa(:,'pelvis_root','head_right',1),-w),dsa(:,'pelvis_root','head_right',1)),h);
ds = circshift(nanmean(GetSegs(ds,1:size(ds,1),lw,nan))',lh);
da = circshift(circ_dist(circshift(dsa(:,'pelvis_root','head_left',1),-w),dsa(:,'pelvis_root','head_left',1)),h);
da = circshift(nanmean(GetSegs(da,1:size(da,1),lw,nan))',lh);
dsv = max([abs(d),abs(ds),abs(da)],[],2);

bs = fet_shake(Trial,'raw');
bs = bs(1:end-(size(bs,1)-dsx.size(1)));

fet = MTADfet('data',[vl(:,'spine_lower')  ,...                         1. [walk] Low pass filtered speed of the lower body
                      vl(:,'spine_middle') ,...                         2. [walk] Low pass filtered speed of the mid body
                      vl(:,'head_front')   ,...                         3. [walk] Low pass filtered speed of the head
                      dsv,...                                           13. smoothed angular speed
                      bs,...                                            14. shake feature
                      dsx(:,'spine_lower',3),...                         4. Height of the lower body
                      dsx(:,'head_front',3)-dsx(:,'spine_lower',3),...   5. Difference in height between head and body
                      dsa(:,'spine_lower','pelvis_root',2) ,...          6. Pitch of the lower back
                      dsa(:,'spine_middle','spine_upper',2),...          7. Pitch of the upper spine
                      ... differential YAW of Joints
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',1),dsa(:,'spine_lower','spine_middle',1))),... 8. [groom] difference in direction of lower back and body
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',1),dsa(:,'spine_lower','spine_upper',1))),... 9. [turn]  difference in direction of lower and upper body
                      abs(circ_dist(dsa(:,'pelvis_root','spine_middle',1),dsa(:,'pelvis_root','spine_upper',1))),...
                      abs(circ_dist(dsa(:,'spine_middle','spine_upper',1),dsa(:,'spine_middle','head_back',1))),...
                      abs(circ_dist(dsa(:,'spine_upper','head_back',1),dsa(:,'spine_upper','head_front',1))),...                      
                      ... differential PITCH of Joints
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',2),dsa(:,'spine_lower','spine_middle',2))),... 
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',2),dsa(:,'spine_lower','spine_upper',2))),... 
                      abs(circ_dist(dsa(:,'pelvis_root','spine_middle',2),dsa(:,'pelvis_root','spine_upper',2))),...
                      abs(circ_dist(dsa(:,'spine_middle','spine_upper',2),dsa(:,'spine_middle','head_back',2))),...
                      abs(circ_dist(dsa(:,'spine_upper','head_back',2),dsa(:,'spine_upper','head_front',2)))],...                      
              'sampleRate', dsx.sampleRate,...
              'syncPeriods',Trial.sync.copy,...
              'syncOrigin', Trial.sync.data(1),...
              'ext','fet',...
              'name','logistic regression feature',...
              'label','fet_lgr',...
              'key','g');
              

