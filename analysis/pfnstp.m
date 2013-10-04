function pfnstp(Session,varargin)
[unit,states,mazeName] = DefaultArgs(varargin,{40,'head','cof'});

%Session = MTASession('jg05-20120310');
%Session = MTATrial(Session,'all')

%% Load MTASession object if Session is type char
if ~isa(Session,'MTASession'),
    Session = MTASession(Session,mazeName);
    Session = MTATrial(Session,'all');
end

%% Load State specific periods
stateLabel = {'theta'};
[stsp,stateLabel] = Session.statePeriods(stateLabel);
stsp = round(stsp./Session.lfpSampleRate.*Session.xyzSampleRate)+1;
%stsp = Session.Bhv.getState(stateLabel).state;
%% Load Units and rescale sampling freq
Session = Session.load_CluRes(Session.xyzSampleRate);
Map = Session.map;
[Res,ind] = SelectPeriods(Session.res,stsp,'d',1,1);
Clu = Session.clu(ind);

%% Select inst XYZ @ spk
stsxyz = SelectPeriods(Session.xyz,stsp, 'c', 1);
stsxyz = reshape(stsxyz,[],size(Session.xyz,2),size(Session.xyz,3));             


%% Select Inst Vel @ spk
rb = Session.Model.rb({'head_back','head_left','head_front','head_right'});
rbcom = sq(Session.com(rb));
rbvel = sqrt(sum(diff(rbcom,1).^2,2)).*Session.xyzSampleRate;
stsrbv = SelectPeriods(rbvel,stsp, 'c', 1);


ang = cell(max(Clu),1);
xyz = cell(max(Clu),1);
xyzvel = cell(max(Clu),1);
spkg = cell(max(Clu),1);
spkgxyz = cell(max(Clu),1);
myspk = cell(max(Clu),1);

for g=1:max(Clu)
    myspk{g} = Res(Clu==g);
    xyz{g} = zeros(length(myspk{g}),size(Session.xyz,2),size(Session.xyz,3));
    xyz{g} = stsxyz(myspk{g},:,:);        
    if length(myspk{g})>10
        isi = diff(myspk{g}/Session.xyzSampleRate*1000);
        isi_tbound = find(isi>100);
        if size(isi_tbound,1)>2,
        isibound = [isi_tbound(1:end),cat(1,isi_tbound(2:end),length(myspk{g}))];
        spkg{g} = zeros(length(isibound),3);
        spkg{g}(:,1) = diff(isibound,1,2);
        for i = 1:length(isibound),
            spkg{g}(i,2) = round(mean([isibound(i,1),isibound(i,2)-1]));
            spkg{g}(i,3) = sum(stsrbv(myspk{g}(isibound(i,1)):myspk{g}(isibound(i,2)-1)))/(length(myspk{g}(isibound(i,1)):myspk{g}(isibound(i,2)))*Session.xyzSampleRate);
            spkg{g}(i,4) = isibound(i,1);
            spkg{g}(i,5) = isibound(i,2)-1;
        end
        spkgxyz{g}   = xyz{g}(spkg{g}(:,2),:,:);
        spkgstart{g} = xyz{g}(spkg{g}(:,4),:,:);
        spkgend{g}   = xyz{g}(spkg{g}(:,5),:,:);
        end
    end
end

% er01-20110721 
% 9 33 44 50 54 60 71

% unit = 33;

tau_limit = 1:2:20;
window_radius = [2:6].*20;

nbins = 100;
pax = -500:(1000/nbins):500;

place_grid = zeros(length(pax),length(pax),length(tau_limit),length(window_radius));
pred_grid = zeros(length(pax),length(pax),length(tau_limit),length(window_radius));
out_grid = zeros(length(pax),length(pax),length(tau_limit),length(window_radius));

for tl = 1:length(tau_limit),
    futime = tau_limit(tl);
    trajdur = round(futime*Session.xyzSampleRate);
    futraj = GetSegs(sq(stsxyz(:,7,:)),myspk{unit}(spkg{unit}(:,4)),trajdur);
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
                place_grid(x,y,tl,wr) = sum(tau==1)/length(tau);
                pred_grid(x,y,tl,wr) = sum(tau>1)/length(tau);
                out_grid(x,y,tl,wr) = sum(tau<0)/length(tau);
            end
        end
        save([Session.spath.analysis Session.filebase '_unit' num2str(unit) '.mat'],'place_grid','pred_grid','out_grid','nbins','pax','window_radius','tau_limit','unit')
    end
end




% 9 33 44
% 50 54 60 71
% $$$ 
% $$$ 
% $$$ 
% $$$ unit = 33;

figure 
t = 1:2:12,
d = 2:6;
for i = 1:length(t),
    for j = 1:length(d),
        %load(['fpf_unit' num2str(unit) '_' num2str(t(i)) 'sec_' num2str(d(j)) 'd.mat']);
        subplotfit(i,j,[length(t),length(d)])
        imagesc(pax,pax,sq(place_grid(:,:,i,j)'))
        %imagesc(pax,pax,sq(pred_grid(:,:,i,j))')
        %imagesc(pax,pax,out_grid)
    end
end

ForAllSubplots('colorbar')





