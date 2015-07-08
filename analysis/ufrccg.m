MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
Trial = MTATrial('jg05-20120317');
Trial.xyz.load(Trial);
Trial.ufr.create(Trial,Trial.xyz);

figure
% $$$ bar(linspace(-round(3*Trial.xyz.sampleRate),round(3*Trial.xyz.sampleRate),round(6*Trial.xyz.sampleRate+1)),mean(reshape(Trial.ufr(round(cell2mat(repmat({Trial.stc{'rear&theta',Trial.xyz.sampleRate}(:,1)},[1,2]))+3*Trial.xyz.sampleRate.*repmat([-1,1],cell2mat({Trial.stc{'rear&theta'}.size(1)}-1),1)),18),[],round(6*Trial.xyz.sampleRate+1)),2));




unit = 29;
state= {'rear&theta','walk&theta'};

yl = [0,0];
sp = [0,0];

figure
for s = 1:numel(state)
sp(s) = subplot2(2,2,1,s); MTAAknnpfs(Trial,unit,state{s}).plot;
end

for s = 1:numel(state)
sp(s) = subplot2(2,2,2,s);
bar(linspace(-3,3,round(6*Trial.xyz.sampleRate)),Filter0(gausswin(21)./sum(gausswin(21)),mean(GetSegs(Trial.ufr(:,unit),round(cell2mat({Trial.stc{state{s},Trial.xyz.sampleRate}(:,1)})-3*Trial.xyz.sampleRate),round(6*Trial.xyz.sampleRate),0),2)));
axis tight
yl(s) = max(ylim)
end
for s = 1:numel(sp),ylim(sp(s),[0,max(yl)]),end



imagesc(GetSegs(Trial.ufr(:,unit),round(cell2mat({Trial.stc{state,Trial.xyz.sampleRate}(:,1)})-3*Trial.xyz.sampleRate),round(6*Trial.xyz.sampleRate),0)');