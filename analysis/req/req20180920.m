% req20180920
%     Status: active
%     Type: Analysis
%     Author: Justin Graboski
%     Final_Forms: NA
%     Project: MjgER2016:interneurons
%     Description: CCG between Interneurons as a function of space
%     Protocol: 


% InterNeuron ccg x head ang

states = {'theta-groom-sit','rear','loc','pause',...
          'groom','sit','R'};

sampleRate = 1250;
Trial = MTATrial.validate('jg05-20120312.cof.all');
xyz = preproc_xyz(Trial,'trb');
xyz.resample(sampleRate);
%ang = create(MTADang,Trial,xyz);
ang = atan2(xyz(:,'head_front',2)-xyz(:,'head_back',2),xyz(:,'head_front',1)-xyz(:,'head_back',1));
stc = Trial.load('stc','msnn_ppsvd_raux');
stcm = stc2mat(stc,xyz,states);
ufr = load(Trial,'ufr',xyz,[],[],0.05,true,'gauss');
lfp = Trial.load('lfp',69);
phz = lfp.phase([30,80]);
phzt = lfp.phase([5,12]);
phz.data = unwrap(phz.data);
phz.resample(xyz);    
phz.data = mod(phz.data+pi,2*pi)-pi;
lfp.resample(xyz);    
vxy = xyz.vel('hcom',[1,2]);
vxy.data(vxy.data<1e-3) = 1e-3;
vxy.data = log10(vxy.data);

unitsInt = select_units(Trial,'int');
ufrInt = Trial.load('ufr',xyz,[],unitsInt,0.04);
int = Trial.load('spk',xyz.sampleRate,'',unitsInt);

binSize = 1;
halfBins = 150;
normalization = 'hz';


center = [200,200];
center = [-200,-200];
center = [0,0];
center = [100,100];
center = [300,0];
center = [-100,-200];
center = [300,-300];

rdist = sqrt(sum(bsxfun(@minus,sq(xyz(:,'hcom',[1,2])),center).^2,2));
%for i = 1:numel(unitsInt)
i = 5;
    iRes = int(unitsInt(i));
    %    for j = i:numel(unitsInt)    
j = 11;    
        jRes = int(unitsInt(j));
sindi = stcm(iRes,1) & ~stcm(iRes,2);
sindj = stcm(jRes,1) & ~stcm(jRes,2);

        
figure,
ny = 5
nx = 5;
angbins = linspace(-pi,pi,nx);
pchbins = linspace(0,2,10);
for d = 1:numel(angbins)-1
    accg1 = [];
    accg2 = [];
    mccg =  [];
% $$$     for e = 1:numel(pchbins)-1
% $$$         grind = sindi &rdist(iRes) < 100;
% $$$         grjnd = sindj & rdist(jRes) < 50;
    
        grind = sindi & (ang(iRes)>angbins(d) & ...
                         ang(iRes)<angbins(d+1)) ...
                & rdist(iRes) < 100;
        grjnd = sindj & (ang(jRes)>angbins(d) & ...
                         ang(jRes)<angbins(d+1)) ...
                & rdist(jRes) < 100;
% $$$         grind = sindi & (ang(iRes,5,7,1)>angbins(d) & ...
% $$$                          ang(iRes,5,7,1)<angbins(d+1)) ...
% $$$                       & (ang(iRes,5,7,2)>pchbins(e) & ...
% $$$                          ang(iRes,5,7,2)<pchbins(e+1)) ...                
% $$$                 & rdist(iRes) < 100;
% $$$         grjnd = sindj & (ang(jRes,5,7,1)>angbins(d) & ...
% $$$                          ang(jRes,5,7,1)<angbins(d+1)) ...
% $$$                 & (ang(jRes,5,7,2)>pchbins(e) & ...
% $$$                    ang(jRes,5,7,2)<pchbins(e+1)) ...                
% $$$                 & rdist(jRes) < 100;
        
% $$$         if sum(grind)==0||sum(grjnd)==0,
% $$$             accg1(:,e) = zeros([2*halfBins+1,1]);
% $$$             accg2(:,e) = zeros([2*halfBins+1,1]);
% $$$             mccg(:,e) =  zeros([2*halfBins+1,1]);
% $$$         else            
            
            [tccg,tbin] = CCG([iRes(grind);jRes(grjnd)],...
                              [ones([sum(grind),1]);2*ones([sum(grjnd),1])],...
                              binSize,halfBins,int.sampleRate,[1,2],normalization);
% $$$             accg1(:,e) = tccg(:,1,1);
% $$$             accg2(:,e) = tccg(:,2,2);        
% $$$             mccg(:,e) =  tccg(:,1,2);        
% $$$         end
% $$$     end
        
        subplot2(ny,nx,1,d);        
        bar(tbin,tccg(:,1,1));
        subplot2(ny,nx,2,d);        
        bar(tbin,tccg(:,1,2));        
        subplot2(ny,nx,3,d);        
        bar(tbin,tccg(:,2,2));        
        subplot2(ny,nx,4,d);                
        plot([phzt(jRes(grjnd));phzt(jRes(grjnd))+2*pi],repmat(circ_dist(phz(jRes(grjnd)), ...
                                                          phz(NearestNeighbour(iRes(grind),jRes(grjnd)))),[2,1]),'.')
        %plot(rdist(jRes(grjnd)),circ_dist(phz(jRes(grjnd)),phz(NearestNeighbour(iRes(grind),jRes(grjnd)))),'.')
        subplot2(ny,nx,5,d);                
        plot([phzt(jRes(grjnd));phzt(jRes(grjnd))+2*pi],repmat(phz(NearestNeighbour(iRes(grind),jRes(grjnd))),[2,1]),'.')
% $$$         subplot2(3,9,1,d);        
% $$$         imagesc(tbin,pchbins,accg1'),axis('xy');
% $$$         subplot2(3,9,2,d);        
% $$$         imagesc(tbin,pchbins,mccg'),axis('xy');
% $$$         subplot2(3,9,3,d);        
% $$$         imagesc(tbin,pchbins,accg2'),axis('xy');        
end
ForAllSubplots('ylim([0,120])');
% $$$ ForAllSubplots('caxis([0,100])');


figure,plot(rdist(jRes(grjnd)),phz(jRes(grjnd)),'.')
figure,plot(rdist(iRes(grind)),phz(iRes(grind)),'.')

figure,plot(phzt(jRes(grjnd)),phz(NearestNeighbour(iRes(grind),jRes(grjnd))),'.')


grind = sindi & rdist(iRes) < 50;
grjnd = sindj & rdist(jRes) < 50;

figure,plot(phzt(jRes(grjnd)),phz(jRes(grjnd)),'.')

figure,plot(rdist(jRes(grjnd)),circ_dist(phz(jRes(grjnd)),phz(NearestNeighbour(iRes(grind),jRes(grjnd)))),'.')
figure,plot(phzt(jRes(grjnd)),circ_dist(phz(jRes(grjnd)),phz(NearestNeighbour(iRes(grind),jRes(grjnd)))),'.')