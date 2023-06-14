

configure_default_args();

MjgER2016_load_data();

trialId = 20;

% RASTER STUFF
Trial = Trials{trialId};    % jg05-20120312.cof.all
unitsSubset = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsPyr    = Trial.spk.get_unit_set(Trial,'pyramidal');
unitsInt =    Trial.spk.get_unit_set(Trial,'interneurons');

% SORT Placefields by spatial bins
pft = pfs_2d_theta(Trial,unitsPyr);
[mxr,mxp]= pft.get_max_rate(unitsPyr);
mxd = sqrt(sum(mxp.^2,2));
mxn = mxp./mxd;
mxa = atan2(mxn(:,2),mxn(:,1));
mxg = nan(size(mxd));
mxg(mxd<=175) = 1;
abins = linspace(-pi,pi,9);
for a = 1:numel(abins)-1
    mxg(mxd>175& mxa>abins(a) & mxa<=abins(a+1)) = a+1;
end    
[mxv,mxi] = sort(mxg,'descend');
unitsPyr = unitsPyr(mxi);
ucolors = hsv(8);
unitsPyrColor = zeros([size(mxg,1),3]);
unitsPyrColor(mxv>1,:) = ucolors(mxv(mxv>1)-1,:);
unitsIntColor = repmat([0,0,0],[numel(unitsInt),1]);

% NEURONQUALITY 
Trial.load('nq');


sampleRate = 30; % Hertz


% LOAD behavioral variables 
xyz = preproc_xyz(Trial,'trb',sampleRate);
fxyz = filter(xyz.copy(),'ButFilter',3,14,'low');
vxy = vel(filter(xyz.copy(),'ButFilter',3,2.5,'low'),{'spine_lower','hcom'},[1,2]);
stc = Trial.stc.copy();


%%%<<< LINEAR 32 lfp spectra -----------------------------------------------------------------------
% LOAD lfp for entire linear shank
channels = 65:96;
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',channels);
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
   
[lys,lfs,lts] = fet_spec(Trial,lfp,[],[],[],specArgsTheta);

% COMPUTE Inter channel theta power correlation
for c1 = 1:32
    for c2 = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c1)),2);
    tp = mean(log10(lys(ind,lfs>6&lfs<12,c2)),2);
    nind = nniz(lp)&nniz(tp);
    crossChanThetaPowCorr(c1,c2) = corr(lp(nind),tp(nind));
    end
end

proxChanThetaPowCorr = [];
tp = mean(log10(cys(ind,cfs>6&cfs<12,2)),2);%+mean(log10(cys(ind,cfs>6&cfs<12,2)),2);
for c = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c)),2);
    nind = nniz(lp)&nniz(tp);
    proxChanThetaPowCorr(c) = corr(lp(nind),tp(nind));
end

%%%>>>


% LOAD lfp for single buszaki shank
Trial.lfp.filename = [Trial.name,'.lfp'];
clfp = Trial.load('lfp',57:64);

% FILTER lfp
fclfp = filter(copy(clfp),'ButFilter',4,30,'low');


elfp = copy(clfp);
elfp.data = [unity(diff(fclfp(:,[1,8]),1,2))];
specArgsTheta = struct('nFFT',2^11,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^10,...
                  'nOverlap',2^10*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[cys,cfs,cts] = fet_spec(Trial,elfp,[],false,[],specArgsTheta);


rcc = [];
tp = mean(log10(cys(ind,cfs>6&cfs<12,2)),2);%+mean(log10(cys(ind,cfs>6&cfs<12,2)),2);
for c = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c)),2);
    nind = nniz(lp)&nniz(tp);
    rcc(c) = corr(lp(nind),tp(nind));
end

% COMPUTE spectra ~ 0.5 second window
elfp = copy(clfp);
elfp.data = [fclfp(:,[1,8])];
specArgsTheta = struct('nFFT',2^10,...
                  'Fs',  lfp.sampleRate,...
                  'WinLength',2^9,...
                  'nOverlap',2^9*0.875,...
                  'NW',3,...
                  'Detrend',[],...
                  'nTapers',[],...
                  'FreqRange',[1,35]);
[ccys,ccfs,ccts,ccphi] = fet_spec(Trial,elfp,[],false,[],specArgsTheta,[],true);



for c1 = 1:32
    for c2 = 1:32
    lp = mean(log10(lys(ind,lfs>6&lfs<12,c1)),2);
    tp = mean(log10(lys(ind,lfs>6&lfs<12,c2)),2);
    nind = nniz(lp)&nniz(tp);
    ccc(c1,c2) = corr(lp(nind),tp(nind));
    end
end

hfig = figure();
imagesc(ccc');
colormap('jet');
caxis([0,1]);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
xlabel(sax(end),'ori          pyr             rad             lm           dg                     ');
ylabel(sax(end),'            dg           lm             rad             pyr          ori');
sax(end+1) = axes(hfig,'Position',[0.15,0.875,0.6,0.025]);
imagesc(rcc);
colormap('jet');
caxis([0,1]);
sax(end).XTickLabel = {};
sax(end).YTickLabel = {};
cax = colorbar(sax(1));
cax.Position(1) = sum(sax(2).Position([1,3]))+0.005;
ylabel(cax,'Theta Power Correlation');
