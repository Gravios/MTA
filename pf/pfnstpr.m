function pfnstpr(Session,varargin)
[states,mazeName] = DefaultArgs(varargin,{'head','cof'});

%Session = MTASession('jg05-20120310');
%Session = MTATrial(Session,'all')

%% Load MTASession object if Session is type char
if ~isa(Session,'MTASession'),
    Session = MTASession(Session,mazeName);
    Session = MTATrial(Session,'all');
end

%% Load State specific periods
[stsp stateLabel] = Session.statePeriods(states);

%% Load Units and rescale sampling freq
[Res Clu Map] = LoadCluRes([Session.spath.nlx Session.name]);
Res = round(Res*Session.lfpSampleRate/Session.sampleRate);
[Res ind] = SelectPeriods(Res,[Session.syncPeriods(1),Session.syncPeriods(end)],'d',1,1);
Clu = Clu(ind);
[myRes ind] = SelectPeriods(Res,stsp,'d',1,1);
myClu = Clu(ind);

%% Select inst XYZ @ spk
stsxyz = SelectPeriods(Session.xyz,round(stsp/Session.lfpSampleRate*Session.xyzSampleRate+1), 'c', 1);
stsxyz = reshape(stsxyz,[],size(Session.xyz,2),size(Session.xyz,3));             

%% Select Inst Vel @ spk
rb = Session.Model.rb({'head_back','head_left','head_front','head_right'});
rbcom = sq(Session.com(rb));
rbvel = sqrt(sum(diff(rbcom,1).^2,2)).*Session.xyzSampleRate;
stsrbv = SelectPeriods(rbvel, round(stsp/Session.lfpSampleRate*Session.xyzSampleRate+1), 'c', 1);

ang = cell(max(Clu),1);
xyz = cell(max(Clu),1);
xyzvel = cell(max(Clu),1);
spkg = cell(max(Clu),1);
spkgxyz = cell(max(Clu),1);
myspk = cell(max(Clu),1);

for g=1:max(Clu)
    myspk{g} = myRes(myClu==g);
    xyz{g} = zeros(length(myspk{g}),size(Session.xyz,2),size(Session.xyz,3));
    xyz{g} = stsxyz(round(myspk{g}/Session.lfpSampleRate*Session.xyzSampleRate+1),:,:);        
    if length(myspk{g})>10
        isi = diff(myspk{g}/Session.lfpSampleRate*1000);
        isi_tbound = find(isi>1000);% miliseconds
        if size(isi_tbound,1)>2,
        isibound = [isi_tbound(1:end),cat(1,isi_tbound(2:end),length(myspk{g}))];
        spkg{g} = zeros(length(isibound),3);
        spkg{g}(:,1) = diff(isibound,1,2);
        for i = 1:length(isibound),
            spkg{g}(i,2) = round(mean([isibound(i,1),isibound(i,2)-1]));
            spkg{g}(i,3) = sum(stsrbv(round(myspk{g}(isibound(i,1))/Session.lfpSampleRate*Session.xyzSampleRate+1):round(myspk{g}(isibound(i,2)-1)/Session.lfpSampleRate*Session.xyzSampleRate+1)))/(length(round(myspk{g}(isibound(i,1))/Session.lfpSampleRate*Session.xyzSampleRate+1):round(myspk{g}(isibound(i,2))/Session.lfpSampleRate*Session.xyzSampleRate+1))*Session.xyzSampleRate);
            spkg{g}(i,4) = isibound(i,1);
            spkg{g}(i,5) = isibound(i,2)-1;
        end
        spkgxyz{g}   = xyz{g}(spkg{g}(:,2),:,:);
        spkgstart{g} = xyz{g}(spkg{g}(:,4),:,:);
        spkgend{g}   = xyz{g}(spkg{g}(:,5),:,:);
        end
    end
end


tau_limit = 1:2:10;
window_radius = [2:4].*20;

nbins = 100;
numBSiterations = 400;
pax = -500:(1000/nbins):500;

place_grid = zeros(length(pax),length(pax),numBSiterations,length(tau_limit),length(window_radius));
pred_grid = zeros(length(pax),length(pax),numBSiterations,length(tau_limit),length(window_radius));
out_grid = zeros(length(pax),length(pax),numBSiterations,length(tau_limit),length(window_radius));

unit = 99;
for bsi = 1:numBSiterations,
    for tl = 1:length(tau_limit),
        futime = tau_limit(tl);
        trajdur = round(futime*Session.xyzSampleRate);
        futraj = GetSegs(sq(circshift(stsxyz(:,7,:),randi([0,round(size(stsxyz,1)/2)]))),round(myspk{unit}(spkg{unit}(:,4))/Session.lfpSampleRate*Session.xyzSampleRate+1),trajdur);
        for wr = 1:length(window_radius),
            tau = -ones(size(futraj,2),1);
            center = [0,0]; 
            for x = 1:length(pax),
                for y = 1:length(pax),
                    center = [pax(x),-pax(y)];
                    tau = -ones(size(futraj,2),1);
                    for i = 1:size(futraj,2),
                        ttau = find(sqrt(sum((sq(futraj(:,i,[1,2]))-repmat(center,trajdur,1)).^2,2))<window_radius(wr),1);
                        if ~isempty(ttau)
                            tau(i) = ttau;
                        else
                            tau(i) = -1;
                        end
                    end
                    place_grid(x,y,bsi,tl,wr) = sum(tau==1)/length(tau);
                    pred_grid(x,y,bsi,tl,wr) = sum(tau>1)/length(tau);
                    out_grid(x,y,bsi,tl,wr) = sum(tau<0)/length(tau);
                end
            end
        end
    end
end

save([Session.spath.analysis Session.filebase '.pfnstpr.' num2str(unit) '.mat'],'place_grid','pred_grid','out_grid','nbins','pax','window_radius','tau_limit','unit')

% $$$ 
% $$$ t = 1:2:20,
% $$$ d = 2:6;
% $$$ figure(301)
% $$$ 
% $$$ set(gcf,'CurrentCharacter','l');
% $$$ unit = 1;
% $$$ while 1,
% $$$     clf
% $$$ 
% $$$ for i = 1:length(t),
% $$$     for j = 1:length(d),
% $$$         %load(['fpf_unit' num2str(unit) '_' num2str(t(i)) 'sec_' num2str(d(j)) 'd.mat']);
% $$$         subplotfit(i,j,[length(t),length(d)])
% $$$         %imagesc(pax,pax,sq(place_grid(:,:,unit,i,j)'))
% $$$         %imagesc(pax,pax,sq(pred_grid(:,:,unit,i,j))')
% $$$         imagesc(pax,pax,sq(out_grid(:,:,unit,i,j)))
% $$$     end
% $$$ end
% $$$ ForAllSubplots('colorbar')
% $$$ 
% $$$     title(num2str(unit))
% $$$     waitforbuttonpress
% $$$     whatkey = get(gcf,'CurrentCharacter');
% $$$     switch double(whatkey)
% $$$       case double('i')
% $$$         unit = input('Enter unit #: ');
% $$$       case double('n')
% $$$         unit = unit+1;
% $$$       case double('p')
% $$$         unit=unit-1;
% $$$       case double('q')
% $$$         return
% $$$     end
% $$$ end
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ 
% $$$ Session = MTASession('er01-20110721');
% $$$ Session = MTATrial(Session,'all');
% $$$ [accg atbin] = autoccg(Session);
% $$$ 
% $$$ Session.Pfs = {};
% $$$ Session = Session.loadPfs;
% $$$ pf_search.stateLabel = 'head';
% $$$ pf_search.mazeName = 'cof';
% $$$ pf_search.trialName = Session.trialName;
% $$$ pf_search.trackingMarker = Session.trackingMarker;
% $$$ pf_search.stateLabel = 'head';
% $$$ pf_search.spk_shuffle = 'n';
% $$$ pf_search.pos_shuffle = 0;
% $$$ pf_search.numBSiterations = 1;
% $$$ pf_search.numZslices = 1;
% $$$ pf_search.nbins = 50;
% $$$ pf_search.smooth = 0.03;
% $$$ 
% $$$ Pfs = Session.getPfs(pf_search);
% $$$ 
% $$$ 
% $$$ figure(302)
% $$$ set(gcf,'CurrentCharacter','l');
% $$$ unit = 1;
% $$$ while 1,
% $$$     clf
% $$$ 
% $$$     subplot2(3,3,1,1)
% $$$     ppf(Pfs.bin1{unit},Pfs.bin2{unit},Pfs.rateMap{unit})
% $$$     subplot2(3,3,1,2)
% $$$     bar(atbin,accg(:,unit));axis tight; grid on
% $$$     subplot2(3,3,1,3)
% $$$     hold on
% $$$     plot(max(sq(max(sq(place_grid(:,:,unit,:,2)),[],1)),[],1),'b')
% $$$     plot(min(sq(min(sq(out_grid(:,:,unit,:,2)),[],1)),[],1),'r')
% $$$     plot(max(sq(max(sq(pred_grid(:,:,unit,:,2)),[],1)),[],1)+max(sq(max(sq(place_grid(:,:,unit,:,2)),[],1)),[],1),'g')
% $$$     ylim([0,1])
% $$$ 
% $$$     subplot2(3,3,2,1)
% $$$     imagesc(pax,pax,sq(out_grid(:,:,unit,2,2))')
% $$$     subplot2(3,3,2,2)
% $$$     imagesc(pax,pax,sq(place_grid(:,:,unit,2,2))')
% $$$     subplot2(3,3,2,3)
% $$$     imagesc(pax,pax,sq(pred_grid(:,:,unit,2,2))')
% $$$     subplot2(3,3,3,1)
% $$$     imagesc(pax,pax,sq(out_grid(:,:,unit,3,2))')
% $$$     subplot2(3,3,3,2)
% $$$     imagesc(pax,pax,sq(place_grid(:,:,unit,3,2))')
% $$$     subplot2(3,3,3,3)
% $$$     imagesc(pax,pax,sq(pred_grid(:,:,unit,3,2))')
% $$$     ForAllSubplots('colorbar')
% $$$ 
% $$$     title(num2str(unit))
% $$$     waitforbuttonpress
% $$$     whatkey = get(gcf,'CurrentCharacter');
% $$$     switch double(whatkey)
% $$$       case double('i')
% $$$         unit = input('Enter unit #: ');
% $$$       case double('n')
% $$$         unit = unit+1;
% $$$       case double('p')
% $$$         unit=unit-1;
% $$$       case double('q')
% $$$         return
% $$$     end
% $$$ end





