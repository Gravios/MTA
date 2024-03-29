
Trial = MTATrial('jg05-20120317');
% $$$ Trial = MTATrial('jg05-20120309');
%Trial = MTATrial('er06-20130612');
%Trial = MTATrial('Ed10-20140812');


xyz = Trial.load('xyz').filter(gtwin(.05,Trial.xyz.sampleRate));

% $$$ [ys,fs,ts] = fet_spec(Trial,xyz.acc(1,3),'mtcsdglong','overwrite',true);
% $$$ [U,S,V] = svd(cov(log10(ys(Trial.stc{'a'},:,1,1))));
% $$$ nfet = zeros([ys.size(1),1]);
% $$$ nfet(nniz(ys)) = log10(ys(nniz(ys),:,1,1))*V(:,1);
% $$$ nfet = MTADxyz('data',nfet,'sampleRate',ys.sampleRate);
% $$$ ufet = reshape([nfet.data';nanmean([nfet.data,circshift(nfet.data,-1)],2)'],1,[])';
% $$$ ufet = MTADxyz('data',reshape([ufet';nanmean([ufet,circshift(ufet,-1)],2)'],1,[])',...
% $$$                'sampleRate',nfet.sampleRate*4);



dsx = xyz.copy;
ufet = dsx;
ufet.resample(30);

dsx.filter(gtwin(.5,dsx.sampleRate));
dsx.resample(ufet);

dsa = Trial.ang.copy;
dsa.create(Trial,dsx);

vel = dsx.vel([],[1,2]);
vel.data(vel.data<0.0001) = 0.0001;
vel.data = log10(vel.data);



% $$$ aper = Trial.sync.copy;
% $$$ aper.resample(xyz.sampleRate);
% $$$ aper.data = aper.data-aper.data(1)+1;
% $$$ aper = aper+[5,-5];
% $$$ 
% $$$ Trial.stc.addState(Trial.spath,...
% $$$                    Trial.filebase,...
% $$$                    aper.data,...
% $$$                    xyz.sampleRate,...
% $$$                    Trial.sync.copy,...
% $$$                    Trial.sync.data(1),...
% $$$                    'gper','a');

dsv = [0;Filter0(ones(1,7)./7,abs(diff(circ_dist(dsa(:,'spine_lower','spine_middle',1),dsa(:,'spine_middle','head_front',1)))))];


wper = resample(Trial.stc{'w'}.cast('TimeSeries'),ufet);
rper = resample(Trial.stc{'r'}.cast('TimeSeries'),ufet);
aper = resample(Trial.stc{'a'}.cast('TimeSeries'),ufet);


lrfet = MTADxyz('data',[vel(:,'spine_lower'),...
                        vel(:,'head_front'),...
                        dsx(:,'spine_lower',3),...
                        dsa(:,'spine_lower','pelvis_root',3),...
                        dsa(:,'spine_middle','spine_upper',2),...
                        dsa(:,'pelvis_root','spine_middle',3),...
                        dsa(:,'spine_lower','pelvis_root',2),...
                        dsa(:,'spine_upper','head_back',3),...
                        dsv,...
                        dsa(:,'spine_lower','head_front',3)],...
                'sampleRate',ufet.sampleRate);
%lrfet.data = [lrfet.data,circshift(lrfet.data,-7),circshift(lrfet.data,7)];


% $$$ lrp = prctile(lrfet(logical(aper.data),:),[5,95])';
% $$$ fbound = zeros([lrfet.size(2),2]);
% $$$ for i=1:lrfet.size(2),
% $$$ fbound(i,:) = [median(lrfet((lrfet(:,i)<lrp(i,1))&aper,i)),median(lrfet((lrfet(:,i)>lrp(i,2))&aper,i))];
% $$$ end
% $$$ lrfet.data = bsxfun(@rdivide,bsxfun(@minus,lrfet.data,fbound(:,1)'),diff(fbound,1,2)')


smat = max(stc2mat(Trial.stc,ufet),[],2);

% smat = [wper(Trial.stc{'a'})+1+rper(Trial.stc{'a'}).*2];
ind = any(smat,2);
B = mnrfit(lrfet(ind,:),smat(ind),'model','nominal');


y = mnrval(B,lrfet.data);

[~,w] = max(y,[],2);

perind =[31600,33700];

figure,
sp(1) = subplot(9,1,[1:4]);
imagesc((1:size(ufet,1))/ufet.sampleRate,1:7,nunity(lrfet(:,1:7))');caxis([0,2]),colormap jet

% $$$ sp(2) = subplot(9,1,5);
% $$$ rper = Trial.stc{'r'};
% $$$ wper = Trial.stc{'w'};
% $$$ rind(1) = find(rper(:,1)>perind(1)-1000,1,'first');
% $$$ rind(2) = find(rper(:,2)<perind(2)+1000,1,'last');
% $$$ wind(1) = find(wper(:,1)>perind(1)-1000,1,'first');
% $$$ wind(2) = find(wper(:,2)<perind(2)+1000,1,'last');
sp(2) = subplot(9,1,[6:9]);
plot((1:size(ufet,1))/ufet.sampleRate,y)
ylim([-.1,1.1])

hold on,Lines(Trial.stc{'w',1}(:),[],'k');
hold on,Lines(Trial.stc{'r',1}(:),[],'r');
hold on,Lines(Trial.stc{'m',1}(:),[],'m');

linkaxes(sp,'x');
xlim(perind/xyz.sampleRate)







% $$$ ind = logical(aper.data);
% $$$ figure,
% $$$ x = abs(circ_dist(dsa(ind,'spine_lower','spine_upper',1),dsa(ind,'spine_middle','head_front',1)));   xl = 0:.04:pi;
% $$$ %x = dsa(ind,'spine_lower','head_front',3);   xl = 20:3:300;
% $$$ y = dsx(ind,'spine_lower',3);   yl = 1:1:70;    
% $$$ %y = dsa(ind,'spine_lower','pelvis_root',2);   yl = .2:.01:1.6;    
% $$$ hist2([x,y],xl,yl),caxis([0,400])
% $$$ 
% $$$ 
% $$$ figure,
% $$$ for i = 1:3,
% $$$     
% $$$     x = circ_dist(dsa(w==i&aper,'spine_lower','spine_upper',1),dsa(w==i&aper,'spine_middle','head_front',1));   xl = 0:.05:pi;
% $$$     %x = dsa(w==i&aper,'spine_lower','pelvis_root',2);   xl = .2:.01:1.6;
% $$$     %x = dsa(w==i&aper,'pelvis_root','spine_middle',2);   xl = -.6:.01:1.4;
% $$$     %x = rhm(w==i&aper.data);   xl = -1.5:.05:1.5;
% $$$     %x = dsa(w==i&aper,'spine_lower','head_front',3);   xl = 20:3:300;
% $$$     %x = ufet(w==i&aper);   xl = 18:.5:65;
% $$$     y = dsa(w==i&aper,'spine_lower','pelvis_root',2);   yl = .2:.01:1.6;    
% $$$     %y = dsx(w==i&aper,'spine_lower',3);                yl = 1:1:70;
% $$$     %y = vel(w==i&aper,'head_front');                    yl = -1.5:.05:2;
% $$$     subplotfit(i,4);
% $$$     hist2([x,y],xl,yl); caxis([0,300])
% $$$ end



Z = [dsa(w==1&aper,'spine_lower','pelvis_root',2),...
     dsa(w==1&aper,'spine_lower','spine_middle',2),...
     dsa(w==1&aper,'spine_lower','head_front',3),...
     dsx(w==1&aper,'spine_lower',3),...
     abs(circ_dist(dsa(w==1&aper,'spine_lower','spine_upper',1),dsa(w==1&aper,'spine_middle','head_front',1)))];


[State, hmm, decode] = gausshmm(nunity(Z),9);


hmm =      hmmtrain(hmm.data.Xtrain,hmm.data.T,hmm);
[decode] = hmmdecode(hmm.data.Xtrain,hmm.data.T,hmm);
State = decode(1).q_star;


ind = w==1&aper;
Stateall = zeros(ufet.size(1),1);
Stateall(ind.data) = State;
figure
plot(Stateall*10+1,'c')
Lines(Trial.stc{'r',ufet.sampleRate}.data(:),[],'r');
Lines(Trial.stc{'w',ufet.sampleRate}.data(:),[],'k');



figure,
for i = 1:9,
    ind = w==1&aper&Stateall==i;
    %ind = w==1&aper;
    %x = dsa(w==1&aper&Stateall==i,'spine_lower','pelvis_root',2);   xl = .2:.01:1.6;
    %x = dsa(ind,'spine_lower','spine_middle',2);          xl = 0:.01:1.2;
    x = abs(circ_dist(dsa(ind,'spine_lower','spine_upper',1),dsa(ind,'spine_middle','head_front',1))); xl = 0:.02:pi;
    %x = ufet(w==i&aper);                                  xl = 18:.5:65;
    %x = dsx(w==i&aper,'spine_lower',3);                   xl = 1:1:70;
    y = dsx(ind,'spine_lower',3);                          yl = 1:1:70;
    %y = dsx(w==1&aper&Stateall==i,'pelvis_root',3);       yl = 15:1:110;
    %y = vel(w==1&aper&Stateall==i,'head_front');          yl = -1.5:.05:2;
    subplotfit(i,9);
    hist2([x,y],xl,yl); caxis([0,200])
end
