% req20180120 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: testing new formulation of rhm
%  Bugs: NA



Trial = MTATrial.validate('Ed01-20140717.cof.all');
Trial = MTATrial.validate('Ed05-20140529.ont.all');
Trial = MTATrial.validate('Ed10-20140817.cof.gnd');
Trial = MTATrial.validate('Ed01-20140709.cof.all');

sessionList = get_session_list('ncp');
Trials = af(@(t) MTATrial.validate(t), sessionList);


xyzProcOpts = 'trb';
xyzProcOpts = '';

eds = linspace(-pi/2,pi/2,30);
theta = linspace(-pi/2,pi/2,30);

mSniff   = nan([numel(Trials),1]);
rotSniff = nan([numel(Trials),1]);
mCohere  = nan([numel(Trials),1]);
rotCohere= nan([numel(Trials),1]);


cpw = nan([numel(Trials),numel(theta),numel(eds)]);
acpw = nan([numel(Trials),numel(theta)]);
scpw = nan([numel(Trials),numel(theta)]);


for tind = 1:numel(Trials),

    Trial = Trials{tind};
    ncpChan = sessionList(tind).ncpChannel;


    rhm = fet_rhm(Trial);
    ncp = fet_ncp(Trial,[],[],ncpChan);
    fet = fet_href_H(Trial,[],[],xyzProcOpts);

    fet.filter('RectFilter');
    fet.data = circshift(fet.data,-1)-fet.data;
    fet.filter('ButFilter',3,[4,20],'bandpass');
    
    fet.data = cat(2,fet.data,rhm.data,ncp.data);

    specArgs = struct('nFFT',2^8,...
                      'Fs',fet.sampleRate,...
                      'WinLength',2^7,...
                      'nOverlap',2^6,...
                      'FreqRange',[1,20]);

    [ysc,fs,ts] = fet_spec(Trial,fet,'mtchglong',false,[],specArgs,'overwrite',true);

    xyz = preproc_xyz(Trial);
    xyz.resample(ysc);
    vxy = xyz.vel('hcom',[1,2]);
    ind = vxy.data>2 & nniz(xyz);
    ang = create(MTADang,Trial,xyz);


    ysn = ysc.copy();
    ysn.data = ysn(:,:,end,end)+eps;
    ysn.filter('RectFilter');
    ysn.data = ysn.data';
    ysn.filter('RectFilter');
    ysn.data = ysn.data';

    [~,sfreqInd] = max(ysn.data,[],2);



    cohereFun = @(y) mean(abs(y(:,1,end)),'omitnan')./mean(sqrt(y(:,1,1).*y(:,end,end)),'omitnan');

    bins = discretize(ang(:,'head_back','head_front',2),eds);

    for t = 1:numel(theta),
        disp(theta(t))
        fet = fet_href_H(Trial,'theta',theta(t));
        fet.filter('RectFilter');
        fet.data = circshift(fet.data,-1)-fet.data;
        fet.filter('ButFilter',3,[4,20],'bandpass');
        
        fet.data = cat(2,fet.data,ncp.data);        

        [ys,fs,ts] = fet_spec(Trial,fet,'mtcsdglong',false,[],specArgs,'overwrite',false);        
        
        ysm = zeros([size(ys,1),2,2]);
        for i = 1:size(ysm,1),
            ysm(i,:,:) = sq(ys(i,sfreqInd(i),[1,end],[1,end]));
        end

        for eid = unique(bins)',
            cpw(tind,t,eid) = cohereFun(ysm(ind&bins==eid&fs(sfreqInd)<14&fs(sfreqInd)>5,:,:));        
        end
        acpw(tind,t) = cohereFun(ysm(ind&fs(sfreqInd)<14&fs(sfreqInd)>5,:,:));
        scpw(tind,t) = log10(median(ysm(ind&fs(sfreqInd)<14&fs(sfreqInd)>5,1,1),'omitnan'));
    end%for t in theta
end%for tind




hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
for tind = 1:numel(Trials),
    subplot2(3,numel(Trials),1,tind);imagesc(theta,eds,sq(cpw(tind,:,:))');
    title(Trials{tind}.filebase);
    subplot2(3,numel(Trials),2,tind);plot(theta,scpw(tind,:)');xlim(theta([1,end]));
    subplot2(3,numel(Trials),3,tind);plot(theta,acpw(tind,:)');xlim(theta([1,end]));
end

for tind = 1:numel(Trials),
    [mSniff(tind),rotSniff(tind)] = max(scpw(tind,:));
    [mCohere(tind),rotCohere(tind)] = max(acpw(tind,:));
end

circ_mean(circ_dist(theta(rotSniff),theta(rotCohere))')

%

tind = 3;




Trial = Trials{tind};
ncpChan = sessionList(tind).ncpChannel;


rhm = fet_rhm(Trial);
ncp = fet_ncp(Trial,[],[],ncpChan);
fet = fet_href_H(Trial,'theta',0);
fet = fet_href_H(Trial,'theta',theta(rotSniff(tind))+0.47);
fet.filter('RectFilter');
fet.data = circshift(fet.data,-1)-fet.data;
fet.filter('ButFilter',3,[4,20],'bandpass');



fet.data = cat(2,fet.data,rhm.data,ncp.data);

specArgs = struct('nFFT',2^8,...
                  'Fs',fet.sampleRate,...
                  'WinLength',2^7,...
                  'nOverlap',2^6,...
                  'FreqRange',[1,20]);



[ysc,fs,ts] = fet_spec(Trial,fet,'mtchglong',true,[],specArgs,'overwrite',true);


hax = gobjects([1,size(ys,3)]);
hfig = figure();
for i = 1:size(ysc,3), 
    hax(end+1)=subplot(size(ysc,3),1,i);
    imagesc(ts,fs,log10(ysc(:,:,i,i))');
    axis('xy');
    colormap('jet');
end;
linkaxes(hax,'xy');

hfig = figure();
for i = 1:size(ys,3), 
    hax(end+1)=subplot(size(ysc,3),1,i);
    if i~=size(ys,3),
        imagesc(ts,fs,ysc(:,:,i,end)');        
        caxis([0.25,1])
    else
        imagesc(ts,fs,log10(ysc(:,:,i,end))');        
    end
    axis('xy');
    colormap('jet');
end;
linkaxes(hax,'xy');
