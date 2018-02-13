


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
units   = cf(@(t)  select_placefields(t),  Trials);

Trials(cell2mat(cf(@isempty,units))) = [];
units(cell2mat(cf(@isempty,units))) = [];

numTrials = numel(Trials);

pft = cf(@(t) pfs_2d_theta(t),  Trials);


pfd = cf(@(t)  MjgER2016_drzfields(t), Trials);
highRateInds = -0.5 < pfd{1}{1}.adata.bins{1} & pfd{1}{1}.adata.bins{1} < 0.5;
% $$$ 
% $$$ nx = 4;
% $$$ ny = 2;
% $$$ t = 18;
% $$$ for u = 1:numel(units{t}),
% $$$ figure(29302302),clf();
% $$$ % PLOT pft
% $$$ subplot2(ny,nx,[1,2],[1,2]);
% $$$ plot(pft{t},units{t}(u),[],true);
% $$$ title(num2str(units{t}(u)));
% $$$ rm = plot(pfd{1}{t},units{t}(u),'isCircular',false);
% $$$ subplot2(ny,nx,[1],3);imagesc(pfd{1}{t}.adata.bins{1},pfd{1}{t}.adata.bins{2},rm');colorbar();axis('xy');
% $$$ subplot2(ny,nx,[2],3);imagesc(pfd{1}{t}.adata.bins{1},pfd{1}{t}.adata.bins{2},bsxfun(@rdivide,rm,max(rm))');colorbar();axis('xy');
% $$$ lrm = LocalMinima2(-rm,0,20);
% $$$ hold('on');
% $$$ scatter(pfd{1}{t}.adata.bins{1}(lrm(1)),pfd{1}{t}.adata.bins{2}(lrm(2)),20,'m','filled');
% $$$ rm = plot(pfd{3}{t},units{t}(u),'isCircular',false);
% $$$ subplot2(ny,nx,[1],4);imagesc(pfd{3}{t}.adata.bins{1},pfd{3}{t}.adata.bins{2},rm');colorbar();axis('xy');
% $$$ subplot2(ny,nx,[2],4);imagesc(pfd{3}{t}.adata.bins{1},pfd{3}{t}.adata.bins{2},bsxfun(@rdivide,rm,max(rm))');colorbar();axis('xy');
% $$$ lrm = LocalMinima2(-rm,0,20);
% $$$ hold('on');
% $$$ scatter(pfd{3}{t}.adata.bins{1}(lrm(1)),pfd{3}{t}.adata.bins{2}(lrm(2)),20,'m','filled');
% $$$ waitforbuttonpress();
% $$$ end

map           = nan([numTrials,max(cellfun('length',units)),2]);
rateMaxInd    = nan([numTrials,max(cellfun('length',units)),numel(pfd{1})]);
rateMaxVal    = nan([numTrials,max(cellfun('length',units)),numel(pfd{1})]);
rateMaxRng    = nan([numTrials,max(cellfun('length',units)),numel(pfd{1}),2]);
rateMaxRngCnt = nan([numTrials,max(cellfun('length',units)),numel(pfd{1})]);


for t = 1:numTrials,
    for u = 1:numel(units{t}),
        for p = 1:numel(pfd{t}),
            rateMap = plot(pfd{t}{p},units{t}(u),'mean',false,[],false,0.95);
            rateMap = mean(rateMap(highRateInds,:),'omitnan');
            [rateMaxVal(t,u,p),rateMaxInd(t,u,p)] = max(rateMap);
            try
                rateMaxRng(t,u,p,:) = cat(4,find([rateMap>[rateMaxVal(t,u,p)/2]]==1,1,'first'),...
                                            find([rateMap>[rateMaxVal(t,u,p)/2]]==1,1,'last'));
                rateMaxRngCnt(t,u,p) = sum(isnan(rateMap(:)));

            end
            map(t,u,:) = cat(3,t,units{t}(u));
        end
    end
end

rateMaxVal    = sq(reshape(rateMaxVal,[],1,numel(pfd{1})));
rateMaxInd    = sq(reshape(rateMaxInd,[],1,numel(pfd{1})));
rateMaxRng    = sq(reshape(rateMaxRng,[],1,numel(pfd{1}),2));
rateMaxRngCnt = sq(reshape(rateMaxRngCnt,[],1,numel(pfd{1})));




binsHPitch  = pfd{1}{1}.adata.bins{2};
binsBPitch = pfd{1}{2}.adata.bins{2};
%binsRHM    = pfd{3}{1}.adata.bins{2};
bins = cf(@(p) p.adata.bins{2}, pfd{1});

rateMax = {};
rateMax{1}  = bins{1}(rateMaxInd(nniz(rateMaxInd),1));
rateMax{2}  = bins{2}(rateMaxInd(nniz(rateMaxInd),2));

rateRngCnt =  rateMaxRngCnt(nniz(rateMaxInd),:,:);


figure();  
for i = 1:2,
    subplot(1,2,i);  hist(rateMax{i}(rateRngCnt(:,i)<20),20);
end

figure();  
ind = rateRngCnt(:,1)<25&rateRngCnt(:,2)<25;
plot(rateMax{2}(ind)+randn([sum(ind),1])/25,...
     rateMax{1}(ind)+randn([sum(ind),1])/25,'.b');
xlim([-pi/2,pi/2]);
ylim([-pi/2,pi/2]);


p = 1;t = 17;figure,for u = 1:numel(units{t}),clf();subplot(1,3,1);plot(pft{t},units{t}(u),'mean',true);subplot(132);plot(pfd{t}{p},units{t}(u),'mean',true,'isCircular',false);subplot(1,3,3);plot(pfd{t}{p+1},units{t}(u),'mean',true,'isCircular',false);waitforbuttonpress();end

subplot(132);  hist(rateMaxHeight(rateRngCnt(:,2)<14),20);
subplot(133);  hist(rateMaxRHM(rateRngCnt(:,3)<15),20);

figure();
subplot(131);  plot(rateMaxPitch+randn([sum(ind),1])/50,...
                    diff(rateMaxRng(nniz(rateMaxRng,1,1),1,:),1,3)+randn(size(rateMaxPitch)),'.');
subplot(132);  plot(rateMaxHeight+randn(size(rateMaxHeight))*10,...
                    diff(rateMaxRng(nniz(rateMaxRng,2,1),2,:),1,3)+randn(size(rateMaxHeight)),'.');
subplot(133);  plot(rateMaxRHM+randn(size(rateMaxRHM))/20,...
                    diff(rateMaxRng(nniz(rateMaxRng,3,1),3,:),1,3)+randn(size(rateMaxRHM)),'.');


rps = 100;
pSc = 1/100;
rSc = 1/20;
hSc = 5;

pSc = 0;
rSc = 0;
hSc = 0;


FigDir = create_directory('/storage/gravio/figures/placefields'); 
hax = gobjects([1,3]);
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hax(1) = subplot(131);  
hold('on');
ind = rateRngCnt(:,1)<12&rateRngCnt(:,2)<6;
out = [];
for i = 1:rps,
    %plot(rateMaxPitch(ind)+randn([sum(ind),1])/25,rateMaxHeight(ind)+randn([sum(ind),1])*10,'.b');
    out = sum(cat(3,out,hist2([rateMaxPitch(ind)+randn([sum(ind),1])*pSc,rateMaxHeight(ind)+randn([sum(ind),1])*hSc],binsPitch,binsHeight)),3);
end
imagesc(binsPitch,binsHeight,out');
axis('tight')
xlabel('pitch (rad)')
ylabel('height (mm)')
title('Peak pitch vs height');


hax(2) = subplot(132);  
hold('on');
ind = rateRngCnt(:,1)<12&rateRngCnt(:,3)<15;
out = [];
for i = 1:rps,
    %plot(rateMaxPitch(ind)+randn([sum(ind),1])/25,rateMaxRHM(ind)+randn([sum(ind),1])/10,'.b');
    out = sum(cat(3,out,hist2([rateMaxPitch(ind)+randn([sum(ind),1])*pSc,rateMaxRHM(ind)+randn([sum(ind),1])*rSc],binsPitch,binsRHM)),3);
end
imagesc(binsPitch,binsRHM,out');
axis('tight')
xlabel('pitch (rad)');
ylabel('rhm 6-12Hz (A.U.)');
title('Peak pitch vs rhm');

hax(3) = subplot(133);
hold('on');
ind = rateRngCnt(:,2)<6&rateRngCnt(:,3)<15;
out = [];
for i = 1:rps,    
    %plot(rateMaxHeight(ind)+randn([sum(ind),1])*10,rateMaxRHM(ind)+randn([sum(ind),1])/10,'.b');
    out = sum(cat(3,out,hist2([rateMaxHeight(ind)+randn([sum(ind),1])*hSc,rateMaxRHM(ind)+randn([sum(ind),1])*rSc],binsHeight,binsRHM)),3);
end
imagesc(binsHeight,binsRHM,out');
axis('tight')
xlabel('height (mm)');
ylabel('rhm 6-12Hz (A.U.)');
title('Peak height vs rhm');
hfig.Position = [0.5,0.5,12,4];
af(@(h) set(h,'Units','centimeters'), hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
FigName = ['pop_drzfields_jpdf'];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));




FigDir = create_directory('/storage/gravio/figures/placefields'); 
hax = gobjects([1,3]);
hfig = figure();
hfig.Units = 'centimeters';
hfig.PaperPositionMode = 'auto';
hax(1) = subplot(131);  
ind = rateRngCnt(:,1)<12&rateRngCnt(:,2)<6;
plot(rateMaxPitch(ind)+randn([sum(ind),1])/25,rateMaxHeight(ind)+randn([sum(ind),1])*10,'.b','MarkerSize',1);
xlabel('pitch (rad)')
ylabel('height (mm)')
title('Peak pitch vs height');

hax(2) = subplot(132);  
ind = rateRngCnt(:,1)<12&rateRngCnt(:,3)<15;
plot(rateMaxPitch(ind)+randn([sum(ind),1])/25,rateMaxRHM(ind)+randn([sum(ind),1])/10,'.b','MarkerSize',1);
xlabel('pitch (rad)');
ylabel('rhm 6-12Hz (A.U.)');
title('Peak pitch vs rhm');

hax(3) = subplot(133);
ind = rateRngCnt(:,2)<6&rateRngCnt(:,3)<15;
plot(rateMaxHeight(ind)+randn([sum(ind),1])*10,rateMaxRHM(ind)+randn([sum(ind),1])/10,'.b','MarkerSize',1);
xlabel('height (mm)');
ylabel('rhm 6-12Hz (A.U.)');
title('Peak height vs rhm');
hfig.Position = [0.5,0.5,12,4];
af(@(h) set(h,'Units','centimeters'), hax);
af(@(h) set(h,'Position',[h.Position(1:2),2,2]), hax);
FigName = ['pop_drzfields'];
print(gcf,'-depsc2',fullfile(FigDir,[FigName,'.eps']));
print(gcf,'-dpng',  fullfile(FigDir,[FigName,'.png']));


