function Stc = bhv_qda(Trial,Stc,varargin)
[train,states,model_filename,display] = DefaultArgs(varargin,{false,[],'MTA_standard_QDA_model.mat',false});


if isempty(states),   states = Trial.stc.list_state_attrib('key'); end
ns = numel(states);

Trial.ang.load(Trial);
Trial.xyz.load(Trial);

dwin= gausswin(61)./sum(gausswin(61));
fwin = gausswin(11)./sum(gausswin(11));
Trial.xyz.filter(fwin);
Trial.ang.load(Trial);

fpv = cat(1,[-2,-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','spine_upper','head_front'},[1,2]))));
fpz = cat(1,[-2,-2],log10(Filter0(dwin,Trial.vel({'spine_lower','head_back'},[3]))));


fet =[];
fet = [fet,fpv,fpz];
fet = [fet,Trial.xyz(:,7,3)-Trial.xyz(:,1,3)];
fet = [fet,Trial.ang(:,3,4,2)];
fet = [fet,Trial.ang(:,5,7,2)];
fet = [fet,fet.^3];
fet = MTADxyz([],[],Filter0(fwin,fet),Trial.xyz.sampleRate);
fet.data(fet.data==0) = nan;
fet_not_nan = ~isnan(fet.data(:,1));


%% Get or Train QDA Model
if train
    Model_Information.description = '';
    Model_Information.StcMode     = Trial.stc.mode;
    Model_Information.StcFilename = Trial.stc.filename;
    Model_Information.SessionFilebase = Trial.filebase;
    Model_Information.state_labels = states;
    keys = Trial.stc.list_state_attrib('key')
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
            spare_keys = 'aqjpxothrymunslf';
            spare_keys = spare_keys(~ismember(spare_keys,cell2mat(keys)));
            Model_Information.state_keys(i) = spare_keys(i);
        end
    end
    save(fullfile(fileparts(mfilename('fullpath')),model_filename),...
         'fet_state','fet_mean_state','cov_state','Model_Information');
    return
end

load(model_filename);
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
             ThreshCross(maxState==i&d_state(:,i)>0,0.5,10),...
             Trial.xyz.sampleRate,...
             Trial.xyz.sync.copy,...
             Trial.xyz.origin,...
             states{i},...
             Model_Information.state_keys(i),...
             'TimePeriods');
end



if display
sts_colors = 'brcmgky';
d_colors   = 'brcmgky';
figure
hold on
for i = 1:ns,
plot(Filter0(dwin,d_state(:,i)),d_colors(i)),
Lines(Trial.stc{states{i},Trial.xyz.sampleRate}(:),[],sts_colors(i));
end
Lines([],0,'k')
end