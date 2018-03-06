% req20180216 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Keywords: Placefields, drz, behavioral transition
%  Description: Computing the mean spike theta phase as a function
%               of drz and behavior onset/offset
%  Bugs: NA


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';
generateFigureParts = false;
marker = 'nose';


if generateFigureParts,
    FigDir = create_directory('/storage/gravio/figures/parts/placefields');
else,
    FigDir = create_directory('/storage/gravio/figures/placefields'); 
end        
        
Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
          cf(@(t)  t.load('nq'),           Trials);

statesTheta = {'loc&theta','lloc&theta','hloc&theta','rear&theta',     ...
               'pause&theta','lpause&theta','hpause&theta',            ...
               'theta-groom-sit'};
states      = {'loc','lloc','hloc','rear','pause','lpause','hpause',...
               'theta-groom-sit'};
numStates = numel(states);

cf(@(t) t.load('nq'), Trials);

t = 1;

Trial = Trials{t};

%ncp = fet_ncp(Trial,[],[],66);

xyz = preproc_xyz(Trial,'trb');                            % LOAD     Marker positions
xyz.filter('RectFilter');                                  % FILTER   marker positions
ang = create(MTADang,Trial,xyz);                           % COMPUTE  Intermarker spherical coordinates

units = select_placefields(Trial);                         % SELECT   subset of hippocampal neurons with high spatial information
pft = pfs_2d_theta(Trial,units);                           % COMPUTE  Expected Neuron firing rate given position in theta state
spk = create( Trial.spk.copy(), Trial, xyz.sampleRate,...  % LOAD     Spike identities and times
             'theta-groom-sit', units, 'deburst');         %                [ theta state; remove bursty spike ( isi < 10ms )]
drz = compute_drz(Trial,units,pft);                        % COMPUTE  Directional Rate Zone
ddz = compute_ddz(Trial,units,pft);                        % COMPUTE  Directional Distance Zone

lfp = Trial.load('lfp',sessionList(t).thetaRef);           % LOAD     Local Field Potential(LFP) of CA1pyr for each electrode shank
lfp.resample(xyz);                                         % RESAMPLE LFP to match xyz
phz = lfp.phase([6,12]);                                   % COMPUTE  LFP theta phase 
%phz = ncp.phase([6,12]);

pch = fet_HB_pitchB(Trial);

vxy = xyz.vel({'spine_lower','hcom'},[1,2]);               %
axy = vxy.copy();                                          %
axy.data = circshift(axy.data,-1)-circshift(axy.data,1);   %
vxy.data(vxy.data<1e-3) = 1e-4;                            % CLIP      lowerbound (<1e-3) of speed to 1e-4
vxy.data = log10(vxy.data);                                % TRANSFORM speed to log10 scale

spkll = create( Trial.spk.copy(), Trial, xyz.sampleRate,...  % LOAD     Spike identities and times
             'lloc&theta', [], 'deburst');         %                [ theta state; remove bursty spike ( isi < 10ms )]

u = 1;

% NCP Stuff
% $$$ figure,
% $$$ rose(phz(spkll(units(u))));
% $$$ [p,th0,r,logZ,k,n] = RayleighTest(phz(spkll.res),spkll.clu);
% $$$ figure,plot(th0,logZ,'.')

% $$$ stsper = Trial.stc{'r'};
% $$$ ststrans = xyz.copy('clear');
% $$$ ststrans.data = nan([size(xyz,1),1]);
% $$$ for period = stsper.data(2:end-1,2)',
% $$$     ststrans.data(period-60:period+60) = linspace(-1,1,121);
% $$$ end

%fet = ststrans;
%fet = pch;
fet = vxy;
binx = linspace([-1,1,10]);                                % SET  bins for x axis drz 
%biny = linspace([-2,pi/2,20]);                             % SET  bins for y axis speed
biny = linspace([-2,2,20]);                                % SET  bins for y axis speed
%biny = linspace([-0.5,0.5,10]);                            % SET  bins for y axis speed
mphzDVSize = [numel(binx),numel(biny)];                    % NOTE final map size
numIter = 1000;                                            % SET  count for bootstrap
figure,
for unit = units
    res = spk(unit);
    res(abs(ddz(res,unit==units))>250) = [];    
    subs = [discretize(drz(res,unit==units),binx),...
            discretize(fet(res,2),biny)];
    %mphz = phz(res,:);    
    mphz = phz(res,spk.map(unit==spk.map(:,1),2));
    mphz(any(isnan(subs),2)) = [];
    subs(any(isnan(subs),2),:) = [];
    if length(subs)>3,
    mphzDV = accumarray(subs,mphz,mphzDVSize,@circ_mean,nan);    
    mphzDVOcc = accumarray(subs,ones([length(subs),1]),mphzDVSize,@sum);
    mphzDVBoot = nan([mphzDVSize,numIter]);
% $$$     for boot = 1:numIter,1,1,1
% $$$         ind = randi(round(0.5*length(subs)),[round(0.5*length(subs)),1]);
% $$$         mphzDVBoot(:,:,boot) = accumarray(subs(ind,:),mphz(ind),mphzDVSize,@circ_mean,nan);
% $$$     end
    clf();
    subplot(1,4,[1,2]);plot(pft,unit,'mean',true,[],false,0.99);
    subplot(143);imagescnan({binx,biny,mphzDVOcc'});axis('xy');
    subplot(144);imagescnan({binx,biny,mphzDV'},[-pi,pi],'circular',true,[0,0,0],1,1,@hsv);axis('xy');
% $$$     subplot(144);imagesc(binx,biny,circ_mean(mphzDVBoot,[],3)');axis('xy');colormap('hsv');colorbar();
    waitforbuttonpress();
    end
end



% $$$ unit = 26;
% $$$ res = spk(unit);
% $$$ figure,
% $$$ %scatter(drz(res,unit==units),ang(res,'spine_middle','spine_upper',2),10,phz(res,spk.map(unit==spk.map(:,1),2)),'filled')
% $$$ scatter(drz(res,unit==units),vxy(res,1),20,phz(res,spk.map(unit==spk.map(:,1),2)),'filled')
% $$$ %scatter(drz(res,unit==units),axy(res,1),20,phz(res,spk.map(unit==spk.map(:,1),2)),'filled')
% $$$ colormap hsv

