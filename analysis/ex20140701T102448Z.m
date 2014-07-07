%ex20140701T102448Z

MTAstartup('mypc','mysd');
Trial = MTATrial('jg05-20120310');
xyz = Trial.xyz.copy;
markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','head_back','head_front'};


fwins = [.25,.5,.75,1,1.25,1.5];


for f = fwins,
    xyz.load(Trial);
    xyz.filter(gtwin(f,xyz.sampleRate));
    vel = xyz.vel(markers,[1,2]);
    for m = 1:numel(markers);
        subplot2(numel(fwins),numel(markers),find(f==fwins),m);
        hist(log10(vel(nniz(vel),m)),100);
        if f == fwins(1), title(markers{m});end
        xlabel(['log10(cm/s) gtwin=' num2str(f)]);
        xlim([-3,3]);
        ylim([0,3e4]);
        set(gca,'YTickLabel',{});
    end
    
end
suptitle('Log10 Vel Distributions');

figure,
f = fwins(1);
vbins = linspace(-2,2,40);
for g = 1:numel(markers);
    for m = 1:numel(markers);
        subplot2(numel(markers),numel(markers),g,m);
        hist2([log10(vel(nniz(vel),m)),log10(vel(nniz(vel),g))],vbins,vbins);
        if  g == 1, title(markers{m});end
        if  m == 1, ylabel(markers{g});end
        %xlabel(['log10(cm/s) gtwin=' num2str(f)]);
%        xlim([-3,3]);
%        ylim([0,3e4]);
%        set(gca,'YTickLabel',{});
    end
end

