

sname = 'jg05-20120317';
marker = 'spine_lower'
stc_mode = 'auto_wbhr';

Trial = MTATrial(sname,'all');
Trial.stc.updateMode(stc_mode);Trial.stc.load;

ang = Trial.ang.copy;
ang.load(Trial);

xyz = Trial.xyz.copy;
xyz.load(Trial);


xyz.filter(gtwin(.5,xyz.sampleRate));
vel = MTADxyz('data',[0;Trial.vel([1],[1,2])],'sampleRate',xyz.sampleRate);


vel.data = log10(vel.data);

ind = ~isinf(vel(:))&~isnan(vel(:))&xyz(:,1,3)~=0;

states = 'rlg';

ind = Trial.stc{'g'};

xbins =  -1:.04:2;
ybins = 20: 15:350;

[N,xbs,ybs] = hist2([vel(ind),xyz(ind,7,3)],xbins,ybins);
N = N./sum(N(:));


[~,hc] = contour(xbs(1:end-1),ybs(1:end-1),N');
set(hc,'LineColor','r')



