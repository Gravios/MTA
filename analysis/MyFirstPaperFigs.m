MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317','all');
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(31)./sum(gausswin(31)));


%% CRAP Figure 1 - Comparison between head and body motion

dmax = 100;dmin =0.1;
dxys = Trial.vel(1:Trial.xyz.size(2),[1,2]);
for i = 1:Trial.xyz.size(2),
    dxysc = clip(dxys(:,i)*10,dmin,dmax);
    dxysc(isnan(dxysc))=dmin;
    [dxyscu(:,i) dxyscux(:,i)] = MakeUniformDistr(dxysc);
end


dxyscu = MTADxyz([],[],dxyscu,Trial.xyz.sampleRate,[],[],Trial.xyz.model);
dxyscux = MTADxyz([],[],dxyscux,Trial.xyz.sampleRate,[],[],Trial.xyz.model);    

tdxyscu = dxyscu.copy;
tdxyscux = dxyscux.copy;
tdxyscu.data = dxyscu(Trial.stc{'t'},:);
tdxyscux.data = dxyscux(Trial.stc{'t'},:);

marpairs = {{{'spine_lower'},{'head_front'}},{{'pelvis_root'},{'head_back'}},{{'spine_upper'},{'head_back'}},{{'head_front'},{'head_back'}}};
nbins = 100;
nticks = 10;
for i = 1:size(marpairs,2),
    figure
    hist2([tdxyscu(:,marpairs{i}{1}) tdxyscu(:,marpairs{i}{2})],100,100)
    caxis([0 100])
    mdxyscux1 = reshape(tdxyscux(1+mod(size(tdxyscux,1),100):end,Trial.xyz.model.gmi(marpairs{i}{1})),[],nbins);
    mdxyscux2 = reshape(tdxyscux(1+mod(size(tdxyscux,1),100):end,Trial.xyz.model.gmi(marpairs{i}{2})),[],nbins);
    for j = 1:nticks
        xticks{j} = sprintf('%2.2f',mdxyscux1(round(size(mdxyscux1,1)/2),j*nticks));
        yticks{j} = sprintf('%2.2f',mdxyscux2(round(size(mdxyscux2,1)/2),j*nticks));
    end
    xl = xlim;
    yl = ylim;
    set(gca,'XTick',linspace(xl(1),xl(2),10),'XTickLabel',xticks)
    set(gca,'YTick',linspace(yl(1),yl(2),10),'YTickLabel',yticks)
end

set(gca,'XTickLabel',sprintf('%2.2f\n', ...
                             mdxyscux1(round(size(mdxyscux1,1)/2),nticks:nticks:nbins)))
set(gca,'YTickLabel',sprintf('%2.2f\n', ...
                             mdxyscux2(round(size(mdxyscux2,1)/2),nticks:nticks:nbins)))
%% End - Figure 1


%% Figure 2 - Marker vs Marker Speed Joint Distribution
MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');

Trial = MTATrial('jg05-20120317','all');
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(31)./sum(gausswin(31)));

states = 'twr';
states = 'wgl';
m1 = 'spine_lower';
m2 = 'head_front';


nbins=100;mab = 2.1;mib = -2;
nbins=50;mab = 2.1;mib = 0.5;

hc = cat(1,mib,mab,mib,mab);
bc = cat(1,mab,mib,mib,mab);


dmax = 100;dmin =0.1;
dxys = Trial.vel(1:Trial.xyz.size(2),[1,2]).*Trial.xyz.sampleRate./10;
dxys = MTADxyz([],[],dxys,Trial.xyz.sampleRate,[],[],Trial.xyz.model);

figure,set(0,'defaulttextinterpreter','none')
for state=states,subplot(3,1,find(state==states));

b = clip(log10(dxys(Trial.stc{state},m1)),mib,mab);
h = clip(log10(dxys(Trial.stc{state},m2)),mib,mab);

hist2([[b;bc],[h;hc]],nbins,nbins);
title([m1 ' vs ' m2 ' during ' Trial.stc{state}.label])
xlabel([m1 ' speed (cm/s)'])
ylabel([m2 ' speed (cm/s)'])
caxis([0,250]);
line([mib,mab],[mib,mab],'color',[1,1,1]);

xticks = get(gca,'XTickLabel');
xticks = mat2cell(xticks,ones(1,size(xticks,1)),size(xticks,2));
xticks = cellfun(@str2num,xticks);
xticks = mat2cell(10.^xticks,ones(1,size(xticks,1)),1);
xticks = cellfun(@sprintf,repmat({'%2.2f'},numel(xticks),1),xticks,'uniformoutput',false);
set(gca,'XTickLabel',xticks);

yticks = get(gca,'YTickLabel');
yticks = mat2cell(yticks,ones(1,size(yticks,1)),size(yticks,2));
yticks = cellfun(@str2num,yticks);
yticks = mat2cell(10.^yticks,ones(1,size(yticks,1)),1);
yticks = cellfun(@sprintf,repmat({'%2.2f'},numel(yticks),1),yticks,'uniformoutput',false);
set(gca,'YTickLabel',yticks);

end


%% End - Figure 2

