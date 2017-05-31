function Stc = bhv_qda(Trial,Stc,varargin)
[train,states,model_filename,display] = DefaultArgs(varargin,{false,[],'MTA_manual_mknsrw_QDA_model.mat',true});


fwin = gtwin(0.75,Trial.xyz.sampleRate);
dwin = gtwin(0.50,Trial.xyz.sampleRate);


xyz = Trial.load('xyz');
xyz.filter(fwin);

ang = Trial.ang.copy; 
ang.create(Trial,xyz);


xyz.addMarker('fhcom',[.7,1,.7],...
              {{'head_back','head_front',[0,0,1]}},...
              xyz.com(xyz.model.rb({'head_left','head_front','head_right'})));

xyz.addMarker('fbcom',[.7,0,.7],...
              {{'spine_lower','pelvis_root',[0,0,1]}},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root'})));


xyz.addMarker('fscom',[.7,0,.7],...
              {{'spine_middle','spine_upper',[0,0,1]}},...
              xyz.com(xyz.model.rb({'spine_middle','spine_upper'})));


fpv = xyz.vel({'fbcom','fscom','fhcom'},[1,2]);
fpv.data(fpv.data<.01) = .01;

fet =[];
fet = [fet,fpv.data];
fet = [fet,xyz(:,'fhcom',3)-xyz(:,'fbcom',3)];
fet = [fet,ang(:,3,4,2)];
fet = [fet,ang(:,5,7,2)];
fet = [fet,ang(:,1,4,3)];
fet = [fet,ang(:,1,3,3)];
fet = [fet,ang(:,2,3,3)];
fet = [fet,fet_turn(Trial)];
fet = [fet,circshift(fet,round(.2*xyz.sampleRate)),circshift(fet,-round(.2*xyz.sampleRate))];

fet = MTADxyz('data',fet,'sampleRate',Trial.xyz.sampleRate);
fet.data(fet.data==0) = nan;
fet_not_nan = prod(~isnan(fet.data),2)==1;



%% Get or Train QDA Model
if train
    
    if isempty(states),   states = Trial.stc.list_state_attrib('label'); end
    ns = numel(states);

    Model_Information.description = '';
    Model_Information.StcMode     = Trial.stc.mode;
    Model_Information.StcFilename = Trial.stc.filename;
    Model_Information.SessionFilebase = Trial.filebase;
    Model_Information.state_labels = states;
    keys = Trial.stc.list_state_attrib('key');
    fet_mean_state = zeros([1,size(fet,2),ns]);
    for i = 1:ns
        fet_state = fet(Trial.stc{states{i}},:);
        fet_state(sum(isnan(fet_state)|...
                      isinf(fet_state),2)>0,:) = [];
        fet_mean_state(1,:,i) = nanmean(fet_state);
        cov_state(:,:,i) = cov(fet_state);
        sti = Trial.stc.gsi(states{i});
        if ~isempty(sti)
            Model_Information.state_keys(i) = keys(sti);
        else
            spare_keys = 'qjpxothrymunslf';
            spare_keys = spare_keys(~ismember(spare_keys,cell2mat(keys)));
            Model_Information.state_keys(i) = spare_keys(i);
        end
    end
    save(fullfile(fileparts(mfilename('fullpath')),model_filename),...
         'fet_state','fet_mean_state','cov_state','Model_Information');
    return
end

load(model_filename);

states = Model_Information.state_labels;
ns = numel(states);

mean_fet_state = repmat(fet_mean_state,[fet.size(1),1,ns]);




%% Transform Features to QDA scores




d_state = zeros(fet.size(1),ns);
for i =  1:ns
    d_state(fet_not_nan,i) = -.5*log(det(cov_state(:,:,i)))...
        -.5*dot(((fet.data(fet_not_nan,:)...
        -mean_fet_state(fet_not_nan,:,i))/cov_state(:,:,i))',(fet.data(fet_not_nan,:)...
        -mean_fet_state(fet_not_nan,:,i))');
end



d_state = Filter0(dwin,d_state);

[~,maxState] = max(d_state,[],2);

for i = 1:ns,
Stc.addState(Trial.spath,...
             Trial.filebase,...
             ThreshCross(maxState==i,0.5,10),...
             Trial.xyz.sampleRate,...
             Trial.xyz.sync.copy,...
             Trial.xyz.origin,...
             states{i},...
             Model_Information.state_keys{i},...
             'TimePeriods');
end


%keyboard
if display
sts_colors = 'brcmgky';
d_colors   = 'brcmgky';
keys = Model_Information.state_keys;
figure
hold on
for i = 1:ns,
plot(Filter0(dwin,d_state(:,i)),d_colors(i)),
%Lines(Trial.stc{keys{i},Trial.xyz.sampleRate}(:),[],sts_colors(i));
end
Lines([],0,'k')
end