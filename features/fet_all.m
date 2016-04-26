function fet = fet_all(Trial,varargin)
%function fet = fet_all(Trial)
%
% exhaustive set of features derived from the raw data
%
%
[newSampleRate,RefTrial] = DefaultArgs(varargin,{20,[]},1);

fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'all_Features','fet_all','a');                  

% Loads preprocessed version of xyz
nm = Trial.xyz.model.N;
[xyz,ss] = preproc_xyz(Trial,'spline_spine');


if ~isempty(RefTrial),
    mfet = fet.copy;
    mfet.data = xyz(:,:,3);
    mfet.map_to_reference_session(Trial,RefTrial);
    xyz.data(:,:,3) = mfet.data;
end              
              
xyz.resample(newSampleRate);

% XYZ filtered 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.5,'low');


% FVELXY Filtered marker speeds in XY plane
fvelxy = xyz.vel([],[1,2]);
fvelxy.filter('ButFilter',3,2.5,'low');
fvelxy.data(fvelxy.data<0)=.1;

facxy = fvelxy.copy;
facxy.data = [diff(facxy.data);zeros([1,fvelxy.size(2)])];

fvelxy.data = log10(fvelxy.data);


% FVELZ Filtered marker speeds in Z axis
fvelz = fxyz.vel([],[3]);
fvelz.filter('ButFilter',3,2.5,'low');
fvelz.data(fvelz.data<0)=.1;

facz = fvelz.copy;
facz.data = [diff(facz.data);zeros([1,fvelz.size(2)])];

fvelz.data = log10(fvelz.data);

% ANG InterMarker Spherical Coordinates
ang = create(MTADang,Trial,xyz);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);


%% Derivatives of features

% ANG (theta) between segments 
d = 2;
cang = ang.copy;
cang.data = [fang(:,1,2,d),...         1
             fang(:,1,3,d),...         2
             fang(:,1,nm+1,d),...      3
             fang(:,1,nm+2,d),...      4
             fang(:,2,3,d),...         5
             fang(:,2,4,d),...         6
             fang(:,3,4,d),...         7
             fang(:,3,nm+2,d),...      8
             fang(:,4,nm+2,d),...      9
             fang(:,4,nm+2,d),...     10
             fang(:,nm+1,nm+2,d),...  11
             fang(:,nm+1,nm+3,d)];   %12
% Add to feature matrix
fet.data =  [fet.data,cang.data];


% DRV of pitch
d = 2;
cang = fang.copy;
cang.data = [ang(:,1,2,d),...         13
             ang(:,1,3,d),...         14
             ang(:,1,nm+1,d),...      15
             ang(:,1,nm+2,d),...      16
             ang(:,2,3,d),...         17
             ang(:,2,4,d),...         18
             ang(:,3,4,d),...         19
             ang(:,3,nm+2,d),...      20
             ang(:,4,nm+2,d),...      21
             ang(:,nm+1,nm+2,d),...   22
             ang(:,nm+1,nm+3,d),...   23
             ang(:,nm+2,nm+3,d)];    %24
cang.filter('ButFilter',3,5,'low');
cang.data = diff(cang.data);
% Add to feature matrix
cang.resample(xyz);
fet.data =  [fet.data,cang.data];



%% DRV of Yaw
d = 1;
cang = ang.copy;
cang.data = [ang(:,1,2,d),...         25
             ang(:,1,3,d),...         26
             ang(:,1,nm+1,d),...      27
             ang(:,1,nm+2,d),...      28
             ang(:,2,3,d),...         29
             ang(:,2,4,d),...         30
             ang(:,3,4,d),...         31
             ang(:,3,nm+2,d),...      32
             ang(:,4,nm+2,d),...      33
             ang(:,4,nm+2,d),...      34
             ang(:,5,nm+2,d),...      35
             ang(:,nm+1,nm+2,d),...   36
             ang(:,nm+1,nm+3,d),...   37
             ang(:,nm+2,nm+3,d)];    %38
cang.filter('ButFilter',3,5,'low');
cang.data = circ_dist(cang.data,circshift(cang.data,1)).*cang.sampleRate;
cang.data = circshift(permute(sum(cang.segs(1:cang.size(1),newSampleRate,0).^2),...
                              [2,3,1]),...
                      [-round(newSampleRate/2),0]);
cang.resample(xyz);
cang.data(cang.data<1e-6) = 1e-6;
cang.data = log10(cang.data);
fet.data =  [fet.data,cang.data];

%% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(fxyz);


%% BFET
mxyz = Trial.load('xyz');
mxyz.filter('ButFilter',3,20);
tsh = round(.15*mxyz.sampleRate);
afet = mxyz.copy;
afet.data = circshift(mxyz.data,-tsh)-circshift(mxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(mxyz(:,4,:)-mxyz(:,1,:),[1,mxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],mxyz.size(2));
afet.resample(fxyz);
zv = afet.copy;
zv.data = log10(abs(afet(:,1))).*sign(afet(:,1));

%% SS
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);
sv = Trial.xyz.copy;
sv.data = sn;
sv.resample(fxyz);


%% AV
sang = [circ_dist(fang(:,1,2,1),fang(:,2,3,1)),...
        circ_dist(fang(:,2,3,1),fang(:,3,4,1)),...
        circ_dist(fang(:,3,4,1),fang(:,4,'hcom',1))];
av = fang.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));


fet.data =  [fet.data,...
             man.data,...        39
             fxyz(:,[1:5],3),... 40,41,42,43,44
             fvelxy(:,{'spine_lower','spine_upper','head_front','bcom','hcom','acom'}),... 45 46 47 48 49 50 
             fvelz(:,{'spine_lower','spine_upper','head_front','bcom','hcom','acom'}),...  51 52 53 54 55 56
             zv.data,... 57
             sv.data,... 58
             av.data...  59
            ];







