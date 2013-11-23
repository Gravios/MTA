Trial = MTATrial('jg05-20120310',{'ang'},'all');
%Trial = MTATrial('jg05-20120310',{'ang',{'lfp',[71,81,91]}},'all');


fuang = Trial.ang(:,4,5,3);
fuang(isnan(fuang)) = mean(fuang(~isnan(fuang)));

%% Filter the trajectories
Trial = Trial.filter();

%% Add head ceter of mass to the xyz field and add markers and
%% sticks to the model
Trial = Trial.addMarker('hcom','[0,0,1]',{{'head_back','hcom',[0,0,1]}},Trial.com(Trial.Model.rb({'head_back','head_left','head_front','head_right'})));

%% Select markers for oscilatory acceleration analysis
rb = Trial.Model.rb({'spine_lower','pelvis_root','spine_upper','hcom'});
sacc = Trial.acc(Trial.Model.gmi(rb.ml));



wacc = WhitenSignal(sacc,[],1);

wind = 2^7;
[ya,fa,ta] = mtchglong(wfuang,2^8,Trial.xyzSampleRate,wind,wind*.875,[],[],[],[.5,50]);
fetSampleRate = 1/diff(ta(1:2));

wdsper = round(Trial.Bhv.getState('walk').state./Trial.xyzSampleRate.*fetSampleRate);

figure,imagesc(1:size(ya,1),fa,(log10(ya)./repmat(log10(max(ya)),size(ya,1),1))'),axis xy
caxis([-8,1.2104])
Lines(wdsper(:,1),[],'k');
Lines(wdsper(:,2),[],'r');

[yw,wind] = SelectPeriods(ya,wdsper,'c',1,1);


wlfp = WhitenSignal(Trial.lfp,[],1);

[yl,fl,tl] = mtchglong(wlfp,2^11,Trial.lfpSampleRate,2^10,512,[],[],[],[1,30]);
spow = mean(ya(:,fa>5&fa<12),2);
tpow = mean(yl(:,fl>5&fl<12),2);
dpow =mean(yl(:,fl<5),2);
dtr=tpow./dpow;

figure,
s1=subplot(311);
imagesc(tg+0.5*tg(1),fg,log10(yg')),axis xy,
s2=subplot(312);
imagesc(ta+0.5*ta(1),fa,log10(ya')),axis xy,caxis([-10,-3.5]),
s3=subplot(313);
imagesc(tl+0.5*tl(1),fl,log10(yl')),axis xy,
linkaxes([s1,s2,s3],'x')



wind =11;
figure
plot(ta+.5*diff(ta(1:2)),Filter0(gausswin(wind)./sum(gausswin(wind)),spow./mean(spow)),'b'),
hold on,
plot(tl+.5*diff(tl(1:2)),Filter0(gausswin(wind)./sum(gausswin(wind)),dtr./mean(dtr)),'m'),
