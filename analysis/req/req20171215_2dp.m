
MjgER2016_load_data();

Trial = Trials{20};
unitSubset = units{20};

% $$$ Trial = MTATrial.validate('jg05-20120310.cof.all');
% $$$ Trial = MTATrial.validate('jg05-20120311.cof.all');
% $$$ Trial = MTATrial.validate('jg05-20120312.cof.all');

pftheta = pfs_2d_theta(Trial,unitSubset);

state = 'pause&theta';
overwrite = false;
sampleRate = 250;



% 2d placefield optimizaiton search
spk = Trial.load('spk',sampleRate,'theta-groom-sit',unitSubset,'deburst');

xyz = Trial.load('xyz','trb');
xyz.resample(sampleRate);
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('hcom', [0.5,1,0.5],{{'head_back','head_front',[0,0,1]}},hcom);
hcom = hcom(:,:,1:2);
xyz.data= xyz.data(:,:,[1,2]);

nx = xyz(:,'head_front',:)-hcom;
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3))); 
ny = nx(:,1,[2,1]);


i = [-300:10:300];
j = [-200:10:200];



% DIAGNOSTIC plot
% $$$ figure();
% $$$ tind = 2110;
% $$$ hold('on');
% $$$ scatter(xyz(tind,{'hcom','head_front','head_right'},1),...
% $$$          xyz(tind,{'hcom','head_front','head_right'},2),...
% $$$         30,[0,0,1;0,1,0;1,0,0],'filled');
% $$$ markerColor = jet(numel(i));
% $$$ for x = 1:numel(i)
% $$$     nxy = sq([nx(tind,1,:)*i(x)+ny(tind,1,:)*j(25)+hcom(tind,1,:)]);
% $$$     scatter(nxy(1),nxy(2),10,markerColor(x,:),'filled');
% $$$ end



% LFP - Local Field Potential
Trial.lfp.filename = [Trial.name,'.lfp'];
try,   lfp = load(Trial,'lfp',sessionList(t).thetaRefGeneral);
catch, lfp = load(Trial,'lfp',sessionList(t).thetaRefGeneral);
end
    
% PHZ - LFP phase within theta band 

phz = lfp.phase([6,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;
lfp.resample(xyz);    


pfd = cell([numel(i),numel(j)]);

for x = 1:numel(i),
    for y = 1:numel(j),
        disp(num2str([x,y]));
        sxyz = bsxfun(@plus,nx*i(x)+ny*j(y),hcom);
        pargs = get_default_args('MjgER2016','MTAApfs','struct');
        pargs.units = unitSubset;
        pargs.states = 'theta-groom-sit-rear';
        pargs.binDims = [40,pi/4,40];
        pargs.SmoothingWeights = [1.5,0.8,1.5];
        pargs.boundaryLimits =  [-500,500;-pi,pi;-500,500];        
        pargs.numIter = 1;
        pargs.overwrite = false;
        pargs.halfsample = false;
        pargs.spk = spk;
        pargs.autoSaveFlag = true;
        pargs.tag = ['ms2dp-x',num2str(i(x)),'y',num2str(j(y))];            
        pargs.xyzp = MTADfet.encapsulate(Trial,[sxyz(:,1,1),phz.data,sxyz(:,1,2)],sampleRate,'','','');
        pargs.compute_pfs = @PlotPFCirc;
        pfsArgs = struct2varargin(pargs);
        pfd{x,y} = MTAApfs(Trial,pfsArgs{:});
    end
end



pff = MTAApfs(Trial,unitSubset,'tag',['ms2dp-x',num2str(i(x)),'y',num2str(j(y))])

bsi = cf(@(p) p.data.si(:,ismember(p.data.clu,unitSubset),:), pfd);
bsit = permute(reshape(cell2mat(bsi),[size(bsi,1),numel(bsi{1}),size(bsi,2),size(bsi,3)]),[1,3,2]);

pmr = cf(@(p) max(p.data.rateMap(:,ismember(p.data.clu,unitSubset),:)), pfd);
pmrt = permute(reshape(cell2mat(pmr),[size(pmr,1),numel(pmr{1}),size(pmr,2),size(pmr,3)]),[1,3,2]);

figure();
for u = 1:size(bsit,3),
    subplot(3,1,1);
    plot(pftheta,unitSubset(u),'mean','text',[],true);
    subplot(3,1,2);    
    imagesc(i,j,pmrt(:,:,u)');
    caxis([0,prctile(reshape(pmrt(:,:,u),[],1),95)])
    axis('xy');
    title(num2str(unitSubset(u)));
    cbx = colorbar();
    cbx.Position = cbx.Position+[0.05,0,0,0];
    subplot(3,1,3);    
    imagesc(i,j,bsit(:,:,u)');
    caxis([0,prctile(reshape(bsit(:,:,u),[],1),95)])    
    axis('xy');
    title(num2str(unitSubset(u)));
    cbx = colorbar();
    cbx.Position = cbx.Position+[0.05,0,0,0];
    waitforbuttonpress();
end

phzOrder = fliplr([5:8,1:4]);

pbins = pfd{1}.adata.bins{2}+double(pfd{1}.adata.bins{2}<0).*2.*pi;
pbins = pbins(phzOrder);
xbins = i([1:16]*4-3);

figure(1);clf();
sp = tight_subplot(8,16,0,0.01);
sp = reshape(sp,16,8)';
for u = 1:size(bsit,3)
    figure(1);
    for b = 1:16,
        rmap = plot(pfd{b*4-3,21},unitSubset(u),[],'text',[],false);
        rmaps = reshape(permute(rmap,[1,3,2]),[],size(rmap,2));
        bsii(b,:) = sum(bsxfun(@rdivide,rmaps,mean(rmaps,'omitnan'))...
                               .*log2(bsxfun(@rdivide,rmaps,mean(rmaps,'omitnan'))),'omitnan')./size(rmaps,1);
        for p = 1:8,
            axes(sp(p,b));
            imagesc(sq(rmap(:,phzOrder(p),:))');
            caxis([0,prctile(rmap(:),99.9)]);
            axis('xy');
            colormap('jet');
            set(gca,'XTicklabel',{});
            set(gca,'YTicklabel',{});
        end
    end
    figure(2);
    imagesc(xbins,pbins,bsii(:,phzOrder)');
    axis('xy');
    waitforbuttonpress();
end



mratesd = cf(@(p,u)  p.maxRate(u),  pfd,repmat({units},size(pfd)));
pareasd = cf(@(p,m)  permute(sum(bsxfun(@gt,p.data.rateMap,m'./2)),[1,3,2]),  pfd,  mratesd);
pd = cell2mat(pareasd);


rd = cf(@(m) permute(m,[2,3,1]),  mratesd);
rd = cell2mat(rd);


FigDir = create_directory('/storage/gravio/figures/placefields_optimization');
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
offset = 12;
hax = gobjects([12,3]);
clf();
for u = 1:12,
% PLOT placefields
    hax(u,1) = subplot2(12,18,u,[1:4]);
    plot(pfd{21,21},units(u+offset));
    if u~=12, clear_axes_ticks(hax(u,1));
    else,     xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('placefield'), end
% PLOT max rate as function of head translation
    hax(u,2) = subplot2(12,18,u,[5:11]);
    imagesc(i,j,sq(rd(:,:,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,2));  
    else,     hax(u,2).YTickLabelMode = 'manual';
              hax(u,2).YTickLabel = {};
              hax(u,2).XTickMode = 'manual';
              hax(u,2).XTick = [-100,-50,0,50,100];
              xlabel('mm');
    end
    if u==1, title('max rate'), end
% PLOT half-max-rate area as function of head translation
    hax(u,3) = subplot2(12,18,u,[12:18]);
    imagesc(i,j,sq(pd(:,:,u+offset))');axis('xy');
    if u~=12, clear_axes_ticks(hax(u,3));
    else,     hax(u,3).YAxisLocation = 'right';
              hax(u,3).XTickMode = 'manual';
              hax(u,3).XTick = [-100,-50,0,50,100];
              xlabel('mm');
              ylabel('mm');
    end
    if u==1, title('half-max-rate area'); end
end
suptitle('max rate and half-max-rate area as function of head translation');
hfig.Position = [0.5,0.5,16,40];
af(@(h)  set(h,'Units','centimeters'),  hax);
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,1));
af(@(h)  set(h,'Position',[h.Position(1:2),2,2]),  hax(:,2:3));
FigName = ['pfs_rate_max_area_min','_2d_',Trial.filebase,'_units-',sprintf('%i-%i',[0,11]+offset+1)];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));





