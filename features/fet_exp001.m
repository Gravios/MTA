function fet = fet_exp001(Trial,varargin)
%function fet = fet_lgr(Trial,sampleRate)
%


if ~isempty(varargin)
sampleRate = varargin{1};
else 
    sampleRate = Trial.xyz.sampleRate;
end

dsx = Trial.load('xyz');

vl = dsx.vel([1,3,7],[1,2]);
vl.filter('ButFilter',3,2,'low');
vl.data(vl.data<0.0001) = 0.0001;
vl.resample(sampleRate);
vl.data = log10(vl.data);

dsx.resample(sampleRate);
dsx.filter('ButFilter',2,4);
dsa = create(Trial.ang.copy,Trial,dsx);


% turning fet smoothed angular speed in the xy plane
dsv = abs(circ_dist(dsa(:,1,4,1),circshift(dsa(:,1,4,1),-1)));


bs = fet_shake(Trial,'raw');
bs = bs(1:end-(size(bs,1)-dsx.size(1)));

fet = MTADfet('data',[vl(:,'spine_lower')  ,...                         1. [walk] Low pass filtered speed of the lower body
                      vl(:,'spine_middle') ,...                         2. [walk] Low pass filtered speed of the mid body
                      vl(:,'head_front')   ,...                         3. [walk] Low pass filtered speed of the head
                      dsv,...                                           4. smoothed angular speed of the spine
                      bs,...                                            5. shake feature
                      dsx(:,'head_front',3)-dsx(:,'spine_lower',3),...  6. Difference in height between head and body
                      ... differential YAW of Joints
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',1),dsa(:,'spine_lower','spine_middle',1))),... 10. [groom] difference in direction of lower back and body
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',1),dsa(:,'spine_lower','spine_upper',1))),...  11. [turn]  difference in direction of lower and upper body
                      abs(circ_dist(dsa(:,'pelvis_root','spine_middle',1),dsa(:,'pelvis_root','spine_upper',1))),... 12.
                      abs(circ_dist(dsa(:,'spine_middle','spine_upper',1),dsa(:,'spine_middle','head_back',1))),...  13.
                      abs(circ_dist(dsa(:,'spine_upper','head_back',1),dsa(:,'spine_upper','head_front',1))),...     14.               
                      ... differential PITCH of Joints
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',2),dsa(:,'spine_lower','spine_middle',2))),... 15.
                      abs(circ_dist(dsa(:,'spine_lower','pelvis_root',2),dsa(:,'spine_lower','spine_upper',2))),...  16.
                      abs(circ_dist(dsa(:,'pelvis_root','spine_middle',2),dsa(:,'pelvis_root','spine_upper',2))),... 17.
                      abs(circ_dist(dsa(:,'spine_middle','spine_upper',2),dsa(:,'spine_middle','head_back',2))),...  18.
                      abs(circ_dist(dsa(:,'spine_upper','head_back',2),dsa(:,'spine_upper','head_front',2)))],...    19.               
              'sampleRate', dsx.sampleRate,...
              'syncPeriods',Trial.sync.copy,...
              'syncOrigin', Trial.sync.data(1),...
              'ext','fet',...
              'name','logistic regression feature',...
              'label','fet_lgr',...
              'key','g');
              

