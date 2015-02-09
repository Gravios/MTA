function Pfs = exp_pfsPerm(Trial,varargin)
[mode,autolabel,niter] = DefaultArgs(varargin,{'compute',false,1000});
states = {{'rear','walk'},{'rear','hswalk'},{'rear','lswalk'},{'hswalk','lswalk'}};


%Trial = MTATrial('jg05-20120317');

if autolabel
    Trial.stc.states = {};
    Trial = labelBhv(Trial);
    Trial = labelAuxBhv(Trial);
    Trial.stc.states{end+1} = theta(Trial);
    Trial.stc.save(1);
else
    Trial.stc.states{Trial.stc.gsi('h')} = Trial.stc{'h'}+[1/Trial.xyz.sampleRate,-1/Trial.xyz.sampleRate];
end




units = select_units(Trial,18,'pyr');
units = units(Trial.nq.SNR(units)>1);
Pfs = {};

switch mode
  case 'compute'
    if matlabpool('size')<4,matlabpool('open',4);end
    parfor s = 1:numel(states),
        for u = units,
            try,MTAApfsPerm(Trial,u,states{s},true,'numIter',niter);end
        end
    end
    if matlabpool('size')==4,matlabpool('close');end
    return
  case 'return'
    units = [];
    for s = 1:numel(states),
        try,Pfs{s} = MTAApfsPerm(Trial,units,states{s},'numIter',niter);end
    end
end


units = [];
aunits=repmat({[]},numel(states),1);;
Pfs=repmat({[]},numel(states),1);
pfs=repmat({[]},numel(states),2);
Tlist = SessionList('pfs');
for i= 1:numel(Tlist),
    Trial = MTATrial(Tlist{i}{1},Tlist{i}{3},Tlist{i}{2});
    Trial.load('nq');
    for s = 1:numel(states),
        pft = MTAApfsPerm(Trial,[],states{s},'numIter',niter);
        units = select_units(Trial,18,'pyr');
        units = units(Trial.nq.SNR(units)>1);
        %pft = MTAApfsPerm(Trial,[],states{1},'numIter',niter);
        units = units(ismember(units,pft.data.clu));
        aunits{s} = [aunits{s},units];
        Pfs{s} = cat(2,Pfs{s},pft.data.rateMap(:,ismember(pft.data.clu,units),:));
        for t = 1:2,
            try,pfs{s,t} = cat(2,pfs{s,t},MTAApfs(Trial,units,states{s}{t}).data.rateMap);end    
        end

    end
end



s =1;
grm = zeros([size(Pfs{s},1),numel(units)]);
for u = 1:numel(units),
srm = sq(Pfs{s}(:,u,:));
for i = 1:size(srm,1),
    try,grm(i,u) = kstest(srm(i,1,:),'alpha',0.001);end
end
end

tpfss = Pfs{s};
Pfs{s}(~repmat(grm,[1,1,size(Pfs{1},3)])) = nan;


s = 4;
zr=[];zw=[];
[z,zind] = max(Pfs{s}(:,:,1));
for u = 1:size(pfs{s,1},2);
zr(u) = pfs{s,1}(zind(u),u,1);
zw(u) = pfs{s,2}(zind(u),u,1);
end

nr=[];nw=[];
[n,nind] = min(Pfs{s}(:,:,1));
for u = 1:size(pfs{s,1},2);
nr(u) = pfs{s,1}(nind(u),u,1);
nw(u) = pfs{s,2}(nind(u),u,1);
end


% $$$ 
% $$$ 
% $$$ figure,plot3(zw,zr,z,'.')
% $$$ figure,plot3(nw,nr,n,'.')
% $$$ 
% $$$ figure,
% $$$ subplot(2,2,1),plot(z,zr,'.')
% $$$ subplot(2,2,2),plot(z,zw,'.')
% $$$ subplot(2,2,3),plot(n,nr,'.')
% $$$ subplot(2,2,4),plot(n,nw,'.')


% $$$ hfig = figure(193939);
% $$$ set(hfig,'paperposition',get(hfig,'paperposition').*[0,0,1,1]+[0,0,4,4])

subplot(121),plot(n,z,'.'),  
%xlim([-55,0]),ylim([0,55]),
xlim([-35,0]),ylim([0,35]),
%title('max z-scores of ratemap differences')
xlabel('max z-scores of the low walk state')
ylabel('max z-scores of the high walk state')

subplot(122),plot(nw,zr,'.'),
%xlim([55,0]),ylim([0,55]),
xlim([0,35]), ylim([0,35])
%title('max rates of placefield ratemaps')
xlabel('max rate of the low walk state (Hz)')
ylabel('max rate of the high walk state (Hz)')




grm = zeros([size(srm,1),numel(units)]);
for u = 1:numel(units),
srm = sq(Pfs{1}(:,u,:));
for i = 1:size(srm,1),
    try,grm(i,u) = kstest(srm(i,1,:));end
end
end


srm(~grm,:) = nan;
u = 50;figure,imagescnan({pft.adata.bins{1},pft.adata.bins{2},reshape(sq(var(srm,[],3)),pft.adata.binSizes')'},[],[],1,[0,0,0]),axis xy,

u = 50;figure,imagesc(pft.adata.bins{1},pft.adata.bins{2},reshape(pfs{1,2}(:,u,:),pft.adata.binSizes')'),axis xy,colorbar







s = 1;
units = [];
for s = 1:numel(states),
    Pfs{s} = MTAApfsPerm(Trial,units,states{s},'numIter',niter);
end

units = Pfs{s}.data.clu;
for t = 1:2,
    pfs{s,t} = MTAApfs(Trial,units,states{s}{t});
end


hfig = figure(74782);
s = 1;
u = units(1);
set(hfig,'position',[2,402,1598,369]);
set(hfig,'paperposition',get(hfig,'paperposition').*[1,1,0,0]+[0,0,12,2.8]);
while u ~= -1,
    mrate = max([max(pfs{1,1}.data.rateMap(:,u==pfs{1,1}.data.clu,:)),max(pfs{1,2}.data.rateMap(:,u==pfs{1,1}.data.clu,:))]);

subplot(131),pfs{s,1}.plot(u,[],1,[0,mrate]);
subplot(132),pfs{s,2}.plot(u,[],1,[0,mrate]);
subplot(133),Pfs{s}.plot(u,'sigks',1);
u = figure_controls(hfig,u,units);
        saveas(hfig,fullfile('/gpfs01/sirota/home/gravio/',...
                             'figures','SFN2014',...
                             ['pfsPerm_' Trial.filebase '-' num2str(u) '.png']),'png');

end



