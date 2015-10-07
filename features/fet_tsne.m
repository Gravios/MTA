function [fet,Nmean,Nstd] = fet_tsne(Trial,varargin)
[newSampleRate,normalized,Nmean,Nstd] = DefaultArgs(varargin,{15,false,[],[]},1);


xyz = Trial.load('xyz');

if isempty(Trial.fet),
    Trial.fet = MTADfet(Trial.spath,...
                        [],...
                        [],...
                        [],...
                        Trial.sync.copy,...
                        Trial.sync.data(1),...
                        []);                  
end


%% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(newSampleRate);


%% XYZ Filtered
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');

%% TRAJ feature
ang = create(MTADang,Trial,xyz);
tsh = 1;
afet = Trial.xyz.copy;
bfet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-tsh)-circshift(fxyz(:,:,[1,2]),tsh);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
[afet.data,bfet.data] = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));
bfet.data = reshape(bfet.data,[],xyz.size(2));
m = MTADxyz('data',circ_dist(afet(:,1),ang(:,1,4,1)),'sampleRate',Trial.xyz.sampleRate);
m.data = circ_dist(circshift(m.data,-5),circshift(m.data,5));
m.data = [diff(m.data);0];
bfet.resample(newSampleRate);

%% XY speed
fvelxy = xyz.vel([],[1,2]);
fvelxy.resample(newSampleRate);
fvelxy.filter('ButFilter',3,2.4,'low');

%% Z speed
fvelz = fxyz.vel([],[3]);
fvelz.resample(newSampleRate);
fvelz.filter('ButFilter',3,2.4,'low');

%% FANG inter marker angles based on filtered xyz
fxyz.resample(newSampleRate);
fang = create(MTADang,Trial,fxyz);



fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'tSNE_Features','fet_tsne','t');                  


fet.data = [fxyz(:,{'spine_lower','spine_upper','head_front'},3),...
            fvelxy(:,{'spine_lower','spine_upper','head_front'}),....
            fvelz(:,'head_back'),...
            man.data,...
            log10(abs(bfet(:,1)+1)),...
            fang(:,'spine_middle','spine_upper',2),...
            fang(:,'spine_upper','head_back',2),...
            fang(:,'head_back','head_front',2),...            
            fang(:,1,4,3).*cos(fang(:,1,4,2)),...
            abs(circ_dist(circshift(fang(:,3,4,2),-1),circshift(fang(:,3,4,2),1))),...
            abs(circ_dist(circshift(fang(:,1,4,1),-1),circshift(fang(:,1,4,1),1))),...
            abs(circ_dist(circshift(fang(:,3,7,1),-1),circshift(fang(:,3,7,1),1)))];
fet.data(isinf(fet(:))) = 0;


if normalized,
    if isempty(Nmean)||isempty(Nstd),
        [~,Nmean,Nstd] = nunity(fet(Trial.stc{'a'},:),@nan);
    end
    [fet.data] = nunity(fet.data,@nan,Nmean,Nstd);
end


