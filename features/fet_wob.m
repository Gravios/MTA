function head_wobble = fet_wob(Trial,varargin)


[sampleRate,mode,windowSize] = DefaultArgs(varargin,{Trial.xyz.sampleRate,'wspectral',1});

% Search for computationally optimal window size 
% given windowSize in seconds Default is 1 second
% swins = 2.^[1:12];
% [~,swi]=min(abs(swins-sampleRate)); 
% windowSize = swins(swi);
%Trial = MTATrial('Ed05-20140528');
Trial = MTATrial('jg05-20120317');

fwin = gtwin(.75,Trial.xyz.sampleRate);
swin = gtwin(.1,Trial.xyz.sampleRate);

xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gausswin(5)./sum(gausswin(5)));

rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
rbb = Trial.xyz.model.rb({'spine_lower','pelvis_root'});

hcom = xyz.com(rb);
Trial.addMarker(xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
Trial.addMarker(xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(fwin,hcom),[1,3,2]));


xyz = Trial.xyz.copy;
xyz.load(Trial);
xyz.filter(gtwin(.1,Trial.xyz.sampleRate));
c = sq(cross(Trial.xyz(:,6,:)-Trial.xyz(:,5,:),Trial.xyz(:,8,:)-Trial.xyz(:,5,:),3));
xyz.filter(gtwin(1,Trial.xyz.sampleRate));
fc = sq(cross(xyz(:,6,:)-xyz(:,5,:),xyz(:,8,:)-xyz(:,5,:),3));

cwob = dot(fc,c,2);
[yso,fso,tso] = mtcsdglong(WhitenSignal(cwob),2^8,Trial.ang.sampleRate,2^7,2^7-1,[],'linear',[],[1,30]);
