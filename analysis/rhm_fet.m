Trial = MTATrial('jg05-20120317');
Trial.load('xyz');
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = Trial.com(rb);
Trial.addMarker(Trial.xyz,'hcom',[.7,0,.7],{{'head_back','head_front',[0,0,1]}},hcom);
Trial.addMarker(Trial.xyz,'fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},permute(Filter0(gausswin(61)./sum(gausswin(61)),hcom),[1,3,2]));

ang = Trial.ang.copy;
ang.create(Trial);
ang = ang(:,5,11,3);
hang(isnan(hang))=nanmean(ang);
%hang = clip(hang,0,20);

wang = WhitenSignal(hang);

%[ya,fa,ta] = mtchglong(diff(Filter0(gausswin(11)./sum(gausswin(11)),hang)),2^9,ang.sampleRate,2^8,2^8*.875,[],[],[],[1,30]);
%[ya,fa,ta] = mtchglong(WhitenSignal(Filter0(gausswin(11)./sum(gausswin(11)),hang)),2^9,ang.sampleRate,2^8,2^8*.875,[],[],[],[1,30]);
[ya,fa,ta] = mtchglong(wang,2^9,ang.sampleRate,2^8,2^8*.875,[],[],[],[1,30]);

figure
imagesc(ta,fa,log10(ya)')
axis xy

