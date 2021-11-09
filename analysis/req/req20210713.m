% req20210713

global MTA_PROJECT_PATH

configure_default_args();

MjgER2016_load_data();

states = {'theta','rear','hloc','hpause','lloc','lpause','groom','sit'};

sampleRate = 250;

phzBins = linspace(0,2*pi,25);

%%%<<< XYHB fine Decoding -------------------------------------------------
tind = [1:28];
cind = 1;
dc = {};
for t = tind;
    Trial = Trials{t};
    tunits = units{t};
    thetaRefChan = sessionList(t).thetaRefGeneral;
    halfSpkWindow = 0.022;
    ufrWindow = 0.012;
    smoothingWeights = [800.^2,800.^2, 1.2.^2, 1.2.^2];    
    %smoothingWeights = [250.^2,250.^2, 0.4.^2, 0.4.^2];
    % THETA
    thetaRefChan = sessionList(t).thetaRefGeneral;
    phzCorr = phzCorrection(t);
    
    xyz = preproc_xyz(Trial,'trb');
    xyz.resample(sampleRate);

    % PLACEFIELDS
    pfs = compute_xyhb_ratemaps( Trial, tunits);
    ds = load(fullfile(MTA_PROJECT_PATH,'analysis','pfsXYHB_mask.mat'));
    mask = ds.mask;

    % UNITS
    spk = Trial.load('spk', sampleRate, '', tunits, '');
    ufr = Trial.load('ufr', xyz,spk,tunits,ufrWindow,'boxcar',true);

    % DECODE
    dc{cind} = decode_ufr_boxcar(Trial, ...
                                 tunits, ...
                                 sampleRate, ...
                                 ufr, ...
                                 pfs, ...
                                 mask,...
                                 halfSpkWindow,...
                                 smoothingWeights,...
                                 'overwrite',true);


    dc{cind}.hRot = headRotation{t};
    
    % THETA PHASE
    dc{cind}.phz = load_theta_phase(Trial,                         ... Trial metadata object
                                    xyz,                           ... Reference object
                                    sessionList(t).thetaRefGeneral,... Theta reference channel
                                    phzCorrection(t));      % Theta phase correction (rad)
    dc{cind}.phz = dc{cind}.phz(dc{cind}.ind,1);
    dc{cind}.iphz = discretize(dc{cind}.phz,phzBins);

    
    dc{cind}.pfs = pfs;
    
    dc{cind}.xyz  = copy(xyz);
    dc{cind}.xyz.data  = xyz.data(dc{cind}.ind,:,:);

    dc{cind}.vxy  = filter(copy(xyz),'ButFilter',4,1.5,'low');;
    dc{cind}.vxy  = dc{cind}.vxy.vel('hcom',[1,2]);
    dc{cind}.vxy.data  = dc{cind}.vxy(dc{cind}.ind,1);
    dc{cind}.lvxy = copy(dc{cind}.vxy);
    dc{cind}.lvxy.data(dc{cind}.lvxy.data<=0) = 0.001;
    dc{cind}.lvxy.data = log10(dc{cind}.lvxy.data);

    dc{cind}.hvec = dc{cind}.xyz(:,'nose',[1,2])-dc{cind}.xyz(:,'hcom',[1,2]);
    dc{cind}.hvec = sq(bsxfun(@rdivide,dc{cind}.hvec,sqrt(sum(dc{cind}.hvec.^2,3))));
    dc{cind}.hvec = cat(3,dc{cind}.hvec,sq(dc{cind}.hvec)*[0,-1;1,0]);
    dc{cind}.hvec = multiprod(dc{cind}.hvec,                           ...
                              [cos(dc{cind}.hRot),-sin(dc{cind}.hRot); ...
                        sin(dc{cind}.hRot), cos(dc{cind}.hRot)],...
                              [2,3],                                   ...
                              [1,2]);

    dc{cind}.tvec = circshift(xyz(:,'hcom',[1,2]),-round(250*0.1)) - ...
                    circshift(xyz(:,'hcom',[1,2]),round(250*0.1));
    dc{cind}.tvec = sq(bsxfun(@rdivide,dc{cind}.tvec,sqrt(sum(dc{cind}.tvec.^2,3))));
    dc{cind}.tvec = cat(3,dc{cind}.tvec,sq(dc{cind}.tvec)*[0,-1;1,0]);
    dc{cind}.tvec = dc{cind}.tvec(dc{cind}.ind,:,:);

    dc{cind}.fet  = fet_HB_pitchB(Trial,dc{cind}.sampleRate);
    dc{cind}.fet.data  = dc{cind}.fet.data(dc{cind}.ind,:);


    dc{cind}.hdist = filter(copy(xyz),'ButFilter',3,5,'low');
    xycoor =    xyz(:,'hcom',[1,2]);
    [~,dc{cind}.hdist.data] = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hdist.data = dc{cind}.hdist.data(dc{cind}.ind);


    % COM 
    dc{cind}.ecom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.com},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.ecom = dc{cind}.ecom{1};
    % SAX 
    dc{cind}.esax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.sax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.esax = dc{cind}.esax{1};
    % MAX 
    dc{cind}.emax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.max},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.emax = dc{cind}.emax{1};
    % LOM 
    dc{cind}.elom = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.lom},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.elom = dc{cind}.elom{1};
    % LAX 
    dc{cind}.elax = cf(@(e,x,h,f)  cat(2,...
                                          sq(multiprod(permute(bsxfun(@minus,...
                                                      e(:,[1,2],:),...
                                                      sq(x(:,'hcom',[1,2]))),...
                                                      [1,2,4,3]),...
                                                      h(:,:,:),2,[2,3])),...
                                          [bsxfun(@minus,f(:,1),e(:,3,:)),bsxfun(@minus,f(:,2),e(:,4,:))]), ...
                          {dc{cind}.lax},{dc{cind}.xyz},{dc{cind}.hvec},{dc{cind}.fet});
    dc{cind}.elax = dc{cind}.elax{1};


    % HVANG
    dc{cind}.hvang = filter(copy(xyz),'ButFilter',4,2,'low');
    xycoor = cat(2,...
                 dc{cind}.hvang(:,'spine_upper',[1,2])-dc{cind}.hvang(:,'bcom',[1,2]),...
                 dc{cind}.hvang(:,'nose',[1,2])-dc{cind}.hvang(:,'hcom',[1,2]));
    dc{cind}.hvang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    % Positive: CCW (Left)     Negative: CW (Right)
    dc{cind}.hvang.data = circ_dist(circshift(dc{cind}.hvang.data(:,2),-10),...
                                    circshift(dc{cind}.hvang.data(:,2),+10));
    dc{cind}.hvang.data = dc{cind}.hvang.data(dc{cind}.ind);

    

% COMPUTE corrected hbang
    dc{cind}.hbang = filter(copy(xyz),'ButFilter',3,30,'low');    
    xycoor = cat(2,...
                 dc{cind}.hbang(:,'spine_upper',[1,2])-dc{cind}.hbang(:,'bcom',[1,2]),...
                 dc{cind}.hbang(:,'nose',[1,2])-dc{cind}.hbang(:,'hcom',[1,2]));
    dc{cind}.hbang.data = cart2pol(xycoor(:,:,1),xycoor(:,:,2));
    dc{cind}.hbang.data = circ_dist(dc{cind}.hbang.data(:,2),dc{cind}.hbang.data(:,1));
    dc{cind}.hbang.data = dc{cind}.hbang.data(dc{cind}.ind);
    dc{cind}.hbang.data = dc{cind}.hbang.data+hbangCorrection(tind(cind));
                       
% GENERATE state collection matrix    
    dc{cind}.stcm = stc2mat(Trial.stc,xyz,states);
    dc{cind}.stcm = dc{cind}.stcm(dc{cind}.ind,:);

% COMPUTE head and body relative velocity vectors
    hvfl = fet_href_HXY(Trial,sampleRate,[],'trb');
    bvfl = fet_bref_BXY(Trial,sampleRate,[],'trb');
    dc{cind}.hvfl = hvfl(dc{cind}.ind,:);
    dc{cind}.bvfl = bvfl(dc{cind}.ind,:);    
    
    cind  = cind +1;
end                     




%%%>>>




%%%<<< CONVERT cell to array -----------------------------------------------------------------------
cind = 1;
dca = dc{cind};
dca.sessionId = tind(cind).*ones(size(dc{cind}.ind));

%dca.xyz = dca.xyz.data;
dca.pos = sq(dca.xyz(:,'hcom',:));
dca.fet = dca.fet.data;
dca.vxy = dca.vxy.data;
dca.lvxy = dca.lvxy.data;
dca.hbang = dca.hbang.data;
dca.hvang = dca.hvang.data;
dca.hdist = dca.hdist.data;

dca.uinc = dc{1}.uinc;
maxUinc = max(cell2array(cf(@(d) size(d.uinc,2),dc)));
dca.uinc = cat(2,dc{1}.uinc,nan(abs([[0,maxUinc]-size(dc{1}.uinc)])));


for cind = 2:numel(dc),
    dca.sessionId = cat(1,dca.sessionId,tind(cind).*ones(size(dc{cind}.ind)));
    dca.ind = cat(1,dca.ind,dc{cind}.ind);
    dca.max = cat(1,dca.max,dc{cind}.max);
    dca.com = cat(1,dca.com,dc{cind}.com);    
    dca.sax = cat(1,dca.sax,dc{cind}.sax);        
    dca.post = cat(1,dca.post,dc{cind}.post);

    
    dca.phz  = cat(1,dca.phz,dc{cind}.phz);    
    dca.iphz = cat(1,dca.iphz,dc{cind}.iphz);    
    
    dca.ucnt = cat(1,dca.ucnt,dc{cind}.ucnt);    
    dca.uinc = cat(1,dca.uinc,cat(2,dc{cind}.uinc,nan(abs([[0,maxUinc]-size(dc{cind}.uinc)]))));        
    %dca.uinc = cat(1,dca.uinc,dc{cind}.uinc);    
    %dca.xyz = cat(1,dca.xyz,dc{cind}.xyz.data);    
    dca.pos = cat(1,dca.pos,sq(dc{cind}.xyz(:,'hcom',:)));
    dca.vxy = cat(1,dca.vxy,dc{cind}.vxy.data);        
    dca.lvxy = cat(1,dca.lvxy,dc{cind}.lvxy.data);            
    
    dca.hvec = cat(1,dca.hvec,dc{cind}.hvec);    
    dca.tvec = cat(1,dca.tvec,dc{cind}.tvec);        
    dca.fet = cat(1,dca.fet,dc{cind}.fet.data);
    
    dca.hvang = cat(1,dca.hvang,dc{cind}.hvang.data);    
    dca.hbang = cat(1,dca.hbang,dc{cind}.hbang.data);        
    dca.hdist = cat(1,dca.hdist,dc{cind}.hdist.data);        
    
    dca.ecom  = cat(1,dca.ecom,dc{cind}.ecom);        
    dca.esax  = cat(1,dca.esax,dc{cind}.esax);        
    dca.emax  = cat(1,dca.emax,dc{cind}.emax);

    dca.elom  = cat(1,dca.elom,dc{cind}.elom);        
    dca.elax  = cat(1,dca.elax,dc{cind}.elax);        
    
    dca.stcm = cat(1,dca.stcm,dc{cind}.stcm);    
    dca.hvfl = cat(1,dca.hvfl,dc{cind}.hvfl);        
    dca.bvfl = cat(1,dca.bvfl,dc{cind}.bvfl);            
    
    dca.hRot  = cat(1,dca.hRot,dc{cind}.hRot);
    dca.pfs   = cat(1,dca.pfs,dc{cind}.pfs);
end
%%%>>>


%%%<<<  Theta(dsc,asc) : hvfl vs ego forward, hbang vs ego lateral
figure();
norm = 'xprob';
%norm = '';
%sind = dca.stcm(:,5)==5 | dca.stcm(:,3)==3;
%sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2 & dca.ucnt>=2;
%dcMode = 'esax';
dcMode = 'elax';
%dcMode = 'ecom';
sind = dca.stcm(:,1)==1 & dca.stcm(:,2)~=2 & ismember(dca.sessionId,[3:5,17:23]) ...
      & dca.hdist < 350;
subplot(221);
    ind = ismember(dca.iphz,[3:10])  & sind;
    hist2([dca.(dcMode)(ind,1),...
           dca.bvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(222);
    ind = ismember(dca.iphz,[16:23]) & sind;
    hist2([dca.(dcMode)(ind,1),...
           dca.bvfl(ind,1)],...
           linspace(-200,200,30),...
           linspace(-10,80,20),norm);
    xlabel('ego forward');
    ylabel('hvfl forward');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(223);
    ind = ismember(dca.iphz,[3:10]) & sind;
    hist2([dca.(dcMode)(ind,2),...
           dca.hbang(ind,1)],...
           linspace(-200,200,30),...
           linspace(-1.5,1.5,20),norm);
    xlabel('ego lateral');
    ylabel('hba rad');
    Lines([],0,'k');
    Lines(0,[],'k');
subplot(224);
    ind = ismember(dca.iphz,[16:23]) & sind;
    hist2([dca.(dcMode)(ind,2),...
           dca.hbang(ind,1)],...
           linspace(-200,200,30),...
           linspace(-1.5,1.5,20),norm);
    xlabel('ego lateral');
    ylabel('hba rad');
    Lines([],0,'k');
    Lines(0,[],'k');

%%%>>>




nbins = 25;
norm = 'xprob';
%norm = '';

figure();
ind = dca.stcm(:,1)==1 & dca.ucnt >= 2  & dca.stcm(:,3)==3 & ismember(dca.sessionId,[3:5,17:23]) & dca.hdist < 350;
hist2([dca.ecom(ind,1),dca.phz(ind)],linspace(-200,300,nbins),phzBins,norm);

figure();
ind = dca.stcm(:,1)==1 & dca.ucnt >= 2  & dca.stcm(:,3)==3 & ismember(dca.sessionId,[3:5,17:23]) & dca.hdist < 350;
hist2([dca.ecom(ind,1),dca.phz(ind)],linspace(-200,300,nbins),phzBins,norm);

figure();
ind =   dca.stcm(:,1) == 1 ...
      & (dca.stcm(:,3) == 3 | dca.stcm(:,5) == 5)...
      & dca.ucnt >= 2  ...
      & ismember(dca.sessionId,[3:5,17:23]) ...
      & dca.hdist < 350;
subplot(121);
hist2([dca.ecom(ind,1),dca.phz(ind)],linspace(-200,300,nbins),phzBins,norm);
colormap('jet')
subplot(122);
hist2([dca.ecom(ind,2),dca.phz(ind)],linspace(-200,200,nbins),phzBins,norm);

figure();
sind = dca.stcm(:,1)==1 & dca.ucnt >= 2 & ~ismember(dca.sessionId,[3:5,17:23]) & dca.hdist < 350;
subplot(121);
    ind = sind & dca.stcm(:,sts)==sts & dca.hbang>0;
    hist2([dca.ecom(ind,2),dca.phz(ind)],linspace(-200,200,nbins),phzBins,norm);
    Lines(0,[],'k');
subplot(122);
    ind = sind & dca.stcm(:,sts)==sts & dca.hbang<0;
    hist2([dca.ecom(ind,2),dca.phz(ind)],linspace(-200,200,nbins),phzBins,norm);
    Lines(0,[],'k');    
colormap('jet')



figure();
sind = dca.stcm(:,1)==1 & dca.ucnt >= 2 & ismember(dca.sessionId,[3:5,17:23]) & dca.hdist < 350;
sind(sind) = logical(randn([sum(sind),1])>0.65);
subplot(121);
    ind = sind & dca.stcm(:,sts)==sts & dca.hbang>0;
    hist2([dca.ecom(ind,2),dca.phz(ind)],linspace(-200,200,nbins),phzBins,norm);
    Lines(0,[],'k');
subplot(122);
    ind = sind & dca.stcm(:,sts)==sts & dca.hbang<0;
    hist2([dca.ecom(ind,2),dca.phz(ind)],linspace(-200,200,nbins),phzBins,norm);
    Lines(0,[],'k');    
colormap('jet')



figure
sts = 5;
ind = dca.stcm(:,1)==1 & dca.stcm(:,sts)==sts & dca.ucnt > 2;
hist2([dca.ecom(ind,3),dca.phz(ind)],linspace(-1.25,1.25,nbins),phzBins,norm);



hfig = open('/storage/share/Projects/BehaviorPlaceCode/decoding/decoding_pos_bhv_part1.fig');
[~,fig,fax,sax] = set_figure_layout(figure(666008),'A4','landscape',[],2.5,2,0,0.2);
hfig.Position(3:4) = [fig.page.width,fig.page.height];
close('666008')
%%%<<< PLOT JPDF ego fwd vs theta phase


% ADJUST subplot coordinates
for sts = 1:6,
    [yind, yOffSet, xind, xOffSet] = deal(7, -1, sts, sts*0.1);
    % CREATE subplot axes
    sax(end+1) = axes('Units','centimeters',                                        ...
                      'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                        fig.page.ypos(yind)+yOffSet,                      ...
                        fig.subplot.width, ...
                        fig.subplot.height],...
                      'FontSize', 8,                                                ...
                      'LineWidth',1);
    hold(sax(end),'on');

    norm = 'xprob';
    nbins = 25;
    ind =   dca.stcm(:,1) == 1 ...
            & dca.stcm(:,sts) == sts ...
            & dca.ucnt >= 2  ...
            & ismember(dca.sessionId,[3:5,17:23]) ...
            & dca.hdist < 350;
    hist2([dca.ecom(ind,1),dca.phz(ind)],linspace(-200,200,nbins),phzBins,norm);
    if sts == 1;
        sax(end).XTick = [-100,0,100];
        sax(end).XTickLabel = {'-10','0','100'};
        xlabel(sax(end),'cm');
    else
        sax(end).XTick = [];
    end
    
    sax(end).YTick = [];
    box(sax(end),'on');
    colormap(sax(end),'jet');
    title(states{sts});

    if sts == 1,
        ylabel('FWD Err')
    end

end
%%%>>>


% ADJUST subplot coordinates

[yind, yOffSet, xind, xOffSet] = deal(7, -1, 7, 0.9);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                    fig.page.ypos(yind)+yOffSet,                                ...
                    fig.subplot.width,                                          ...
                    fig.subplot.height],                                        ...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

norm = 'xprob';
nbins = 25;
ind =   dca.stcm(:,1) == 1 ...
        & (dca.stcm(:,7) ~= 7 | dca.stcm(:,8) ~= 8) ...
        & dca.ucnt >= 2  ...
        & ismember(dca.sessionId,[3:5,17:23]) ...
        & dca.hdist < 350;
hist2([dca.ecom(ind,3),dca.phz(ind)],linspace(-1.25,1.25,nbins),phzBins,norm);
sax(end).XTick = [-0.5,0,0.5];
sax(end).YTick = [];
box(sax(end),'on');
colormap(sax(end),'jet');
Lines(0,[],'m');
title(sax(end),{'theta','HP Err'});

%%%>>>


% ADJUST subplot coordinates

[yind, yOffSet, xind, xOffSet] = deal(7, -1, 8, 1.0);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,                      ...
                    fig.page.ypos(yind)+yOffSet,                      ...
                    fig.subplot.width, ...
                    fig.subplot.height],...
                  'FontSize', 8,                                                ...
                  'LineWidth',1);
hold(sax(end),'on');

norm = 'xprob';
nbins = 25;
ind =   dca.stcm(:,1) == 1 ...
        & (dca.stcm(:,7) ~= 7 | dca.stcm(:,8) ~= 8) ...
        & dca.ucnt >= 2  ...
        & ismember(dca.sessionId,[3:5,17:23]) ...
        & dca.hdist < 350;
hist2([dca.ecom(ind,4),dca.phz(ind)],linspace(-1.25,1.25,nbins),phzBins,norm);
sax(end).XTick = [];
sax(end).YTick = [];
box(sax(end),'on');
colormap(sax(end),'jet');
Lines(0,[],'m');
title(sax(end),{'theta','BP Err'});

% ADJUST subplot coordinates        
[yind, yOffSet, xind, xOffSet] = deal(7, -1, 9, 1.1);
% CREATE subplot axes
sax(end+1) = axes('Units','centimeters',                        ...
                  'Position',[fig.page.xpos(xind)+xOffSet,      ...
                    fig.page.ypos(yind)+yOffSet,      ...
                    fig.subplot.width/1.5,              ...
                    fig.subplot.height],              ...
                  'FontSize', 8,                                ...
                  'LineWidth',1);
hold(sax(end),'on');
% ADD phase "bar"
plot(-cos(circ_ang2rad(-60:420)),-60:420,'-k','LineWidth',2);
sax(end).YAxisLocation = 'right';
sax(end).Visible = 'off';                                    
sax(end).YTick = [];
sax(end).YTickLabels = [];
sax(end).XTick = [-0.5,0,0.5];
ylim(sax(end),[-60,420]);
xlim(sax(end),[-1.5,1.5]);

text( 0.25, 360, '2\pi', 'HorizontalAlignment','center',     ...
      'FontSize', 8, 'VerticalAlignment',  'middle');
text(-0.25, 180, '\pi', 'HorizontalAlignment','center',     ...
     'FontSize', 8, 'VerticalAlignment',  'middle');                    
text( 0.25,   0,   '0', 'HorizontalAlignment','center',     ...
      'FontSize', 8, 'VerticalAlignment',  'middle');                    
%%%>>>
savefig(hfig,'/storage/share/Projects/BehaviorPlaceCode/decoding/decoding_pos_bhv_part2.fig');

% $$$ 
% $$$ 
% $$$ norm ='yprob';
% $$$ mthd = 'esax';
% $$$ %mthd = 'ecom';
% $$$ figure,
% $$$ for sts = 1:6,
% $$$     sind = dca.stcm(:,1)==1 & dca.stcm(:,sts)==sts & ismember(dca.sessionId,[3:5,17:23]);
% $$$     subplot2(6,4,sts,1);    
% $$$     ind = sind & dca.iphz==10;
% $$$     hist2([dca.ucnt(ind),...
% $$$            dca.(mthd)(ind,3)],...
% $$$           linspace(0,30,31),...
% $$$           linspace(-1.2,1.2,30),norm);
% $$$     title(states{sts});
% $$$ 
% $$$ 
% $$$     subplot2(6,4,sts,2);
% $$$     hist2([dca.ucnt(ind),...
% $$$            dca.(mthd)(ind,4)],...
% $$$           linspace(0,30,31),...
% $$$           linspace(-1.2,1.2,30),norm);
% $$$     title(states{sts});
% $$$ 
% $$$     subplot2(6,4,sts,3);    
% $$$     ind = sind& dca.iphz==18;
% $$$     hist2([dca.ucnt(ind),...
% $$$            dca.(mthd)(ind,3)],...
% $$$           linspace(0,30,31),...
% $$$           linspace(-1.2,1.2,30),norm);
% $$$     title(states{sts});
% $$$ 
% $$$ 
% $$$     subplot2(6,4,sts,4);
% $$$     hist2([dca.ucnt(ind),...
% $$$            dca.(mthd)(ind,4)],...
% $$$           linspace(0,30,31),...
% $$$           linspace(-1.2,1.2,30),norm);
% $$$     title(states{sts});
% $$$ 
% $$$ end
% $$$ 
% $$$ ForAllSubplots('caxis([0,0.1])');