% req20170829 
% spine spline svd

Trial = MTATrial.validate('jg05-20120317.cof.all');
Trial = MTATrial.validate('jg05-20120309.cof.all');
varargin = {};

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('xyzProcOpts',    {'SPLINE_SPINE_HEAD_EQD'},                                    ...
                 'stcMode',        'msnn_ppsvd'                                                  ...
);
[xyzProcOpts,stcMode] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

stc      = load(Trial,'stc',stcMode);
[xyzs,ss] = preproc_xyz(Trial,xyzProcOpts);
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
angs = create(MTADang,Trial,xyzs);



% Tranlational movements relative to body
shft = 3;
tvec = zeros([size(ss,[1,2]),2]);
for m = 1:size(ss,2),
    tvec(:,m,:) = circshift(ss(:,m,[1,2]),-shft)-circshift(ss(:,m,[1,2]),shft);
end

unvec = [];
rotationAngles = deg2rad([0,90]);
mvec = xyz(:,'spine_upper',[1,2])-xyz(:,'spine_lower',[1,2]);
for theta = rotationAngles,
    rotMat = repmat(permute([cos(theta),-sin(theta);sin(theta),cos(theta)],[3,1,2]),[size(mvec,1),1,1]);
    unvec(:,end+1,:) = bsxfun(@rdivide,multiprod(mvec,rotMat,[2,3],[2,3]),sqrt(sum(mvec.^2,3)));
end

nind = nniz(tvec);
dwalkFetRot = zeros([size(nind,1),numel(rotationAngles),size(ss,2)]);
for t = rotationAngles;
    for m = 1:size(ss,2),
        dwalkFetRot(nind,t==rotationAngles,m) = dot(tvec(nind,m,:),unvec(nind,t==rotationAngles,:),3);
    end
end
dwalkFetRot = permute(dwalkFetRot,[1,3,2]);

ssr = ss.copy();
ssr.data(:,:,[1,2]) = dwalkFetRot;


[~,S,V] = svd(reshape(ssr(stc{'a'},:,:),[],prod(size(ssr,[2,3]))),0);


figure();
subplot(6,1,1);
plot(diag(S));
for i = 1:5,
    subplot(6,1,i+1);
    sV = reshape(V(:,i+10),[],3);
    plot(sV);
end



spineSegmentLength = MTADxyz('data',sqrt(sum(diff(ss.data,1,2).^2,3)),'sampleRate',xyzs.sampleRate);
spineLengthXY = MTADxyz('data',sqrt(sum(diff(ss.data(:,:,[1,2]),1,2).^2,3)),'sampleRate',xyz.sampleRate);

sper = stc{'a-s'};
sper.cast('TimeSeries');
sper.data = logical(sper.data);


figure();histogram(ang(stc{'a-s'},1,2,3),linspace(40,120,100))
figure();histogram(sum(spineLength(:,1:77),2),linspace(40,120,100))


% ss(t,m,d)
nind = find(nniz(xyz));
ssn = ss.copy();
ssn.data = zeros([size(ss,1),size(ss,2)-6,size(ss,3)]);
markerInd = zeros([size(ss,1),3])
for t = nind'
    nzind = [true,spineSegmentLength(t,:)]>1e-6;
    zind = find(~nzind)-1;
    nzind([zind]) = ~nzind([zind]);    
    ssn.data(t,:,:) = ss(t,nzind,:);
    markerInd(t,:) = zind-[1:numel(zind)].*2;
end

xyzn = Trial.load('xyz');

numMarkers = 5;
baseInd = 100/(numMarkers-1);


for m = 1:numMarkers-2,
    medianMarkerIndOffset = baseInd*m-median(markerInd(nind,m));
    for t = nind'
        xyzn.data(t,m,:) = ssn(t,markerInd(nind,m)+medianMarkerIndOffset,:)
    end
end

spineSegmentLengthNew = MTADxyz('data',sqrt(sum(diff(ssn.data,1,2).^2,3)),'sampleRate',xyz.sampleRate);



figure,imagesc(spineSegmentLengthNew(:,:)'),


[xyzn,ssn] = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQI_SMOOTH');



Trial = MTATrial('jg05-20120317');
Trial = MTATrial('er01-20110719');
Trial = MTATrial('Ed01-20140707');
Trial.load('stc','hand_labeled_rev2_jg');
Trial.load('stc','hand_labeled_rev3_jg');

Trials = af(@(t) MTATrial.validate(t), get_session_list('BHV_S4H5'));
Trials = af(@(t) MTATrial.validate(t), get_session_list('hand_labeled'));
f = cf(@(t) list_files(t.name,'.sehs.'), Trials);

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,1.5);
ang = create(MTADang,Trial,xyz);

figure();
hold('on');
ind = Trial.stc{'s&a'};plot(xyz(ind,3,3),ang(ind,1,4,2),'.b')
ind = Trial.stc{'p&a'};plot(xyz(ind,3,3),ang(ind,1,4,2),'.c')
ind = Trial.stc{'m&a'};plot(xyz(ind,3,3),ang(ind,1,4,2),'.m')




er = 5
jg = 9
xyz_er = Trials{er}.load('xyz','sehs'); 
xyz_jg = Trials{jg}.load('xyz','sehs'); 
xyzo_er = Trials{er}.load('xyz'); 
xyzo_jg = Trials{jg}.load('xyz'); 

%xyzo = Trials{11}.load('xyz'); 
%xyzc = Trials{11}.load('xyz','sehs'); 

stc_er = Trials{er}.load('stc','msnn_ppsvd');
stc_jg = Trials{jg}.load('stc','msnn_ppsvd');

ang_er = create(MTADang,Trials{er},xyz_er);
ang_jg = create(MTADang,Trials{jg},xyz_jg);
ango_er = create(MTADang,Trials{er},xyzo_er);
ango_jg = create(MTADang,Trials{jg},xyzo_jg);


eds = linspace(20,180,100);

state = 'w+p+n+r';
figure,
subplot(211);hold('on');
bar(eds,histc(xyz_er(stc_er{state},4,3),eds),'histc');
hax = bar(eds,histc(xyzo_er(stc_er{state},4,3),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.5;
hax.EdgeAlpha = 0.5;
subplot(212);hold('on');
bar(eds,histc(xyz_jg(stc_jg{state},4,3),eds),'histc');
hax = bar(eds,histc(xyzo_jg(stc_jg{state},4,3),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.5;
hax.EdgeAlpha = 0.5;
drawnow;


eds = linspace(0,180,100);
state = 'w+p+n+r';
figure,
subplot(211);hold('on');
bar(eds,histc(ang_er(stc_er{state},2,4,3),eds),'histc');
hax = bar(eds,histc(ango_er(stc_er{state},2,4,3),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.5;
hax.EdgeAlpha = 0.5;
subplot(212);hold('on');
bar(eds,histc(ang_jg(stc_jg{state},2,4,3),eds),'histc');
hax = bar(eds,histc(ango_jg(stc_jg{state},2,4,3),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.5;
hax.EdgeAlpha = 0.5;
drawnow;


eds = linspace(-0.8,pi/2,100);
state = 'w+p+n+r';
figure,
subplot(211);hold('on');
bar(eds,histc(ang_er(stc_er{state},2,4,2),eds),'histc');
hax = bar(eds,histc(ango_er(stc_er{state},2,4,2),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.5;
hax.EdgeAlpha = 0.5;
subplot(212);hold('on');
bar(eds,histc(ang_jg(stc_jg{state},2,4,2),eds),'histc');
hax = bar(eds,histc(ango_jg(stc_jg{state},2,4,2),eds),'histc');
hax.FaceColor = 'r';
hax.EdgeColor = 'r';
hax.FaceAlpha = 0.5;
hax.EdgeAlpha = 0.5;
drawnow;

figure();hold('on');
plot(xyzo(:,4,3));
plot(xyzc(:,4,3));