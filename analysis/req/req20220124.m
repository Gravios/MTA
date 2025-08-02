
MjgER2016_load_data();	  
Trial = Trials{23};		  
xyz = preproc_xyz(Trial,'trb');

fxyz = filter(copy(xyz),'ButFilter',4,2.4,'low');

hfxyz = filter(copy(xyz),'ButFilter',4,20,'low');

vxy = vel(fxyz,{'spine_lower','hcom'},[1,2]);


fvxy = vel(fxyz,{'spine_lower','hcom'},[1,2]);
fvxy.data(fvxy.data<1e-4) = 1e-4;
fvxy.data = log10(fvxy.data);

figure();
hold('on');
plot(fvxy(:,1))
plot(lvxy(:,1))

stc = Trial.load('stc','hand_labeled_rev3_jg');


figure();
v = 2;
subplot(211);hold('on');
ind = stc{'w'};
set(histogram(fvxy(ind,2)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'p'};
set(histogram(fvxy(ind,2)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'n'};
set(histogram(fvxy(ind,2)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'s'};
set(histogram(lvxy(ind,2)),'EdgeColor','none','FaceAlpha',0.5);

figure();
sts = 'wrnpgs';
for s = 1:numel(sts)
subplot2(6,2,s,1);
    ind = stc{sts(s)};
    hist2(fvxy(ind,1:2),linspace(-4,2,50),linspace(-4,2,50));
    line([-4,2],[-4,2],'Color','w');
subplot2(6,2,s,2);
    ind = stc{sts(s)};
    hist2(lvxy(ind,1:2),linspace(-4,2,50),linspace(-4,2,50));
    line([-4,2],[-4,2],'Color','w');
end


fang = create(MTADang,Trial,fxyz);

figure,plot(  fang(:,'spine_lower','pelvis_root',3) ...
            + fang(:,'pelvis_root','spine_middle',3) ...
            + fang(:,'spine_middle','spine_upper',3) ...
            + fang(:,'spine_upper','hcom',3))

sleng = copy(xyz);
sleng.data =   fang(:,'spine_lower','pelvis_root',3) ...
            + fang(:,'pelvis_root','spine_middle',3) ...
            + fang(:,'spine_middle','spine_upper',3) ...
            + fang(:,'spine_upper','hcom',3);

[xyz,ss] = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');


ang = create(MTADang,Trial,hfxyz);

figure,
subplot(211);
plot([1:size(sv,1)]./sv.sampleRate,sv.data)
subplot(212);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


hxy = fet_href_HXY(Trial,[],[],'trb',20);




figure,
subplot(211);
hold('on');
plot([1:size(hxy,1)]./hxy.sampleRate,hxy(:,[1]))
plot([1:size(hxy,1)]./hxy.sampleRate,hxy(:,[2]))
plot([1:size(hxy,1)]./hxy.sampleRate,hxy(:,[2])+hxy(:,[1]))
%plot([1:size(hxy,1)]./hxy.sampleRate,hxy(:,[3]))
%plot([1:size(hxy,1)]./hxy.sampleRate,circ_dist(ang(:,'hcom','nose',1),circshift(ang(:,'hcom','nose',1),-1))*120)
subplot(212);
    hold('on');
    plotSTC(Trial.stc,1,'text',{'theta','sit','groom','lpause','lloc','hpause','hloc','rear'},'kymbbggr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

x = hxy(:,[2]);
x = newx(:,2);
x(~nniz(x))=0;

[ys,fs,ts] = mtchglong(WhitenSignal(x),2^9,hxy.sampleRate,2^8,2^8.*0.875,[],[],[],[0,20]);

[ys,fs,ts] = mtchglong(x,2^9,hxy.sampleRate,2^8,2^8.*0.875,[],[],[],[0,20]);
ts =ts+diff(ts(1:2))/2;


xts = [1:size(hang.roll)]./xyz.sampleRate;
figure,
subplot(211);hold('on');
%imagesc(ts,fs,log10(bsxfun(@rdivide,ys,mean(ys(:,fs<1),2)))')
imagesc(ts,fs,log10(ys)')
% $$$ plot(ts,mean(log10(bsxfun(@rdivide,ys(:,fs>3&fs<6),mean(ys(:,fs<1),2))),2)')
% $$$ plot(ts,mean(log10(bsxfun(@rdivide,ys(:,fs>10&fs<15),mean(ys(:,fs<1),2))),2)')
% $$$ plot(ts,mean(log10(bsxfun(@rdivide,ys(:,fs>15),mean(ys(:,fs<1),2))),2)')
% $$$ plot(xts,circ_dist(hroll,circshift(hroll,-10)));
% $$$ plot(xts,circ_dist(hpitch,circshift(hpitch,-10)));
% $$$ plot(xts,circ_dist(hang.direction,circshift(hang.direction,-10)))
% $$$ legend({'roll','pitch','yaw'});
axis('xy');
colormap('jet');
subplot(212); 
    hold('on');
    plotSTC(stc,1,'text',{'theta','sit','groom','pause','turn','walk','rear'},'kymcgbr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');


%570*xyz.sampleRate
tstart = 68332;
%563*xyz.sampleRate
tstart = 67493;
figure
ind = tstart-200:tstart;
plot3(ind(1:end-1),...
      diff(hpitch(ind)),...
      circ_dist(hang.direction(ind(1:end-1)),hang.direction(ind(2:end))));



figure,
hist2([hang.roll(1:end-1),diff(hpitch)],linspace(-pi,pi,100),linspace(-0.06,0.06,100));
caxis([0,1000])

dpitch = circ_dist(hpitch,circshift(hpitch,-1));
dyaw = circ_dist(hyaw,circshift(hyaw,-1));
droll = circ_dist(hroll,circshift(hroll,-1));

py = multiprod([dyaw,dpitch],[cos(pi/4),-sin(pi/4);sin(pi/4),cos(pi/4)],2,[1,2]);
py = multiprod([dyaw,dpitch],[cos(pi/3),-sin(pi/3);sin(pi/3),cos(pi/3)],2,[1,2]);

figure,
subplot(211);
hold('on');
plot(xts,py(:,1))
plot(xts,py(:,2))
plot(xts,dyaw)
colormap('jet');
subplot(212); 
    hold('on');
    plotSTC(stc,1,'text',{'theta','sit','groom','pause','turn','walk','rear'},'kymcgbr');
    ylim([1,9]);
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,2)/5+1,'r');
    plot([1:size(vxy)]./vxy.sampleRate,vxy(:,1)/5+1,'g');
    xlabel('Time (s)');
linkaxes(findobj(gcf(),'Type','Axes'),'x');

norm = 'xprob'
norm = '';
figure,
subplot(311);
    ind = cast(stc{'m'},'TimeSeries');
    ind = logical(ind.data);
    hist2([hang.roll(ind)-0.4,dyaw(ind)],linspace(-1,1,50),linspace(-0.06,0.06,50),norm);
subplot(312);
    ind = cast(stc{'n'},'TimeSeries');
    ind = logical(ind.data);
    hist2([hang.roll(ind)-0.4,dyaw(ind)],linspace(-1,1,50),linspace(-0.06,0.06,50),norm);
subplot(313);
    ind = cast(stc{'w'},'TimeSeries');
    ind = logical(ind.data);
    hist2([hang.roll(ind)-0.4,dyaw(ind)],linspace(-1,1,50),linspace(-0.06,0.06,50),norm);

hang.roll(ind)-0.4,dyaw(ind)
    
    

hang  = transform_origin(Trial,hfxyz,'head_back','head_front',{'head_left','head_right'});
hroll = hang.roll;
hroll(nniz(hroll)) = ButFilter(hang.roll(nniz(hroll)),4,[1,20]./(xyz.sampleRate/2),'bandpass');
hpitch = hang.pitch;
hpitch(nniz(hpitch)) = ButFilter(hang.pitch(nniz(hroll)),4,[1,20]./(xyz.sampleRate/2),'bandpass');
hyaw = hang.direction;
hyaw(nniz(hyaw)) = ButFilter(hang.direction(nniz(hyaw)),4,[20]./(xyz.sampleRate/2),'low');


%log10(bsxfun(@rdivide,ys,max(ys,[],2)))'





figure();
hold('on');
ind = stc{'w'};
set(histogram(sleng(ind,1)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'p'};
set(histogram(sleng(ind,1)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'n'};
set(histogram(sleng(ind,1)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'s'};
set(histogram(sleng(ind,1)),'EdgeColor','none','FaceAlpha',0.5);
ind = stc{'m'};
set(histogram(sleng(ind,1)),'EdgeColor','none','FaceAlpha',0.5);
