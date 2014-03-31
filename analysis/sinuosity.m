function sinos = sinuosity(Data,varargin);
[markers,dims,windowSize,filterWindow] = DefaultArgs(varargin,{[1:8],[1,2,3],121,[]});
% Diagnostic Args
%Trial = MTATrial('jg05-20120310','all');
%fwin = gausswin(11)./sum(gausswin(11));
%mind = [1:8];
%dim = 3;
%wins = 241;
if Data.isempty, error('MTAData:sinuosity data field is empty');end
if isempty(filterWindow),filterWindow = gausswin(9)./sum(gausswin(9));end

xyz = Data.copy;
xyz.filter(filterWindow);

xsegs = GetSegs(sq(xyz(:,markers,dims)),1:xyz.size(1),windowSize,0);
xsegs = reshape(xsegs,windowSize,xyz.size(1),numel(markers),numel(dims));

vsegs = sum(sqrt(sum(diff(xsegs).^2,4)));
dsegs = sqrt(sum((xsegs(end,:,:,:)-xsegs(1,:,:,:)).^2,4));

tsinos = sq(vsegs./dsegs)-1;

sinos = xyz;
sinos.filename = [];
sinos.clear;
sinos.data = tsinos;
if sinos.size(1)==1,sinos.data = sinos.data';end

% $$$ 
% $$$ vel = [zeros(1,numel(mind));Trial.vel(mind,dim)];
% $$$ vel = MTADxyz('data',vel,'sampleRate',Trial.xyz.sampleRate);
% $$$ vel.filter(gausswin(61)./sum(gausswin(61)));
% $$$ 
% $$$ figure,plot(log10(sins(:,1)));
% $$$ hold on,plot(log10(vel(:,1)*5),'r')
% $$$ 
% $$$ 
% $$$ figure
% $$$ e1 = linspace(-2,2,60);
% $$$ e2 = linspace(-1.8,1.8,60);
% $$$ m = 1;
% $$$ subplot(331),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0,m),sinos(vel(:,m)~=0,m)]),-2,2),e1,e2)
% $$$ m = 3;
% $$$ subplot(332),
% $$$ hist2(clip(log10([sinos(vel(:,m)~=0,1),sinos(vel(:,m)~=0,m)]),-2,2),e2,e2)
% $$$ m = 7;
% $$$ subplot(333),
% $$$ hist2(clip(log10([sinos(vel(:,m)~=0,1),sinos(vel(:,m)~=0,m)]),-2,2),e2,e2)
% $$$ m = 2;
% $$$ subplot(334),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0,m),sinos(vel(:,m)~=0,1)]),-2,2),e1,e2)
% $$$ m = 4;
% $$$ subplot(335),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0,m),sinos(vel(:,m)~=0,1)]),-2,2),e1,e2)
% $$$ m = 7;
% $$$ subplot(336),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0,m),sinos(vel(:,m)~=0,1)]),-2,2),e1,e2)
% $$$ m = 3;
% $$$ subplot(337),
% $$$ hist2(clip(log10([vel(vel(:,1)~=0,1),sinos(vel(:,1)~=0,m)]),-2,2),e1,e2)
% $$$ m = 4;
% $$$ subplot(338),
% $$$ hist2(clip(log10([vel(vel(:,1)~=0,1),sinos(vel(:,1)~=0,m)]),-2,2),e1,e2)
% $$$ m = 7;
% $$$ subplot(339),
% $$$ hist2(clip(log10([vel(vel(:,1)~=0,1),sinos(vel(:,1)~=0,m)]),-2,2),e1,e2)



% $$$ 
% $$$ rp = Trial.stc{'r'}.copy;rp.cast('TimeSeries');
% $$$ %rp.data = ~rp.data;
% $$$ rp.data = ones(rp.size);
% $$$ e1 = linspace(-2,2,60);
% $$$ e2 = linspace(-1.5,3,60);
% $$$ 
% $$$ figure
% $$$ m = 1;
% $$$ subplot(331),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0&rp(:),m),sinos(vel(:,m)~=0&rp(:),m)]),-2,5),e1,e2)
% $$$ m = 3;
% $$$ subplot(332),
% $$$ hist2(clip(log10([sinos(vel(:,m)~=0&rp(:),1),sinos(vel(:,m)~=0&rp(:),m)]),-2,5),e2,e2)
% $$$ m = 7;
% $$$ subplot(333),
% $$$ hist2(clip(log10([sinos(vel(:,m)~=0&rp(:),1),sinos(vel(:,m)~=0&rp(:),m)]),-2,5),e2,e2)
% $$$ m = 2;
% $$$ subplot(334),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0&rp(:),m),sinos(vel(:,m)~=0&rp(:),1)]),-2,5),e1,e2)
% $$$ m = 4;
% $$$ subplot(335),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0&rp(:),m),sinos(vel(:,m)~=0&rp(:),1)]),-2,5),e1,e2)
% $$$ m = 7;
% $$$ subplot(336),
% $$$ hist2(clip(log10([vel(vel(:,m)~=0&rp(:),m),sinos(vel(:,m)~=0&rp(:),1)]),-2,5),e1,e2)
% $$$ m = 3;
% $$$ subplot(337),
% $$$ hist2(clip(log10([vel(vel(:,1)~=0&rp(:),1),sinos(vel(:,1)~=0&rp(:),m)]),-2,5),e1,e2)
% $$$ m = 4;
% $$$ subplot(338),
% $$$ hist2(clip(log10([vel(vel(:,1)~=0&rp(:),1),sinos(vel(:,1)~=0&rp(:),m)]),-2,5),e1,e2)
% $$$ m = 7;
% $$$ subplot(339),
% $$$ hist2(clip(log10([vel(vel(:,1)~=0&rp(:),1),sinos(vel(:,1)~=0&rp(:),m)]),-2,5),e1,e2)
% $$$ 
% $$$ 
% $$$ 
% $$$ figure,
% $$$ m = 1;
% $$$ hist2(clip(log10([Trial.xyz(vel(:,1)~=0&rp(:),7,3),sinos(vel(:,1)~=0&rp(:),m)]),-2,3),100,100)
% $$$ hist2([clip(log10(Trial.xyz(vel(:,1)~=0,7,3)),1.5,3),clip(log10(sinos(vel(:,1)~=0,m)),-1.8,2.8)],100,100)
% $$$ 
% $$$ 
% $$$ m = 1;
% $$$ figure,
% $$$ subplot(311),bar(e2,histc(log10(sinos(Trial.stc{'r'},m)),e2),'histc');
% $$$ subplot(312),bar(e2,histc(log10(sinos(Trial.stc{'w'},m)),e2),'histc');
% $$$ subplot(313),bar(e2,histc(clip(log10(sinos(sinos(:,3)~=0,3)),e2(1),e2(end)),e2),'histc');
% $$$ ForAllSubplots('ylim([0,30000])')
% $$$ ForAllSubplots('Lines([],10000,''k'')')
