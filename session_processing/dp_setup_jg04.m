%% MTASession Setup - jg04 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite = false;
nlx = true;
report = false;
trb = false;
S = get_session_list('jg04_CA1');
S = get_session_list('jg04_CA3');
S = get_session_list('Ed10');
i = 2;

if overwrite, QuickSessionSetup(S(i)); end
s = MTASession.validate([S(i).sessionName,'.',S(i).mazeName,'.',S(i).trialName]); 

% IF you forgot to include a vsk upon construction
% s.xyz.model =  MTAModel(fullfile(s.spath, [s.name '-' s.maze.name '.vsk']),'-vsk');
% s.xyz.save



% PLOT Session creation report
hfig = figure(gen_figure_id),
sp(1) = subplot2(6,4,1:2,1:2);  pZ(s);        % PLOT timeseries for z axis
sp(2) = subplot2(6,4,1:2,3:4);  pXY(s);       % PLOT position within xy plane
sp(3) = subplot2(6,4,3:4,1:4);  pSE(s,hfig);  % PLOT simple check for rigidbody errors

% CORRECT marker swaps 
headMarkers = regexpi(s.xyz.model.ml,'^head_[A-Za-z]*','match');
headMarkers = cellfun(@(x) cellstr(x),headMarkers(~cellfun('isempty',headMarkers)));
% automatic
ERCOR_fillgaps_RidgidBody(s,'EMGM','BEST_SWAP_PERMUTATION',[],headMarkers);
% manual
ERCOR_fillgaps_RidgidBody(s,'MANUAL','BEST_SWAP_PERMUTATION',[],headMarkers);

% RECONSTRUCT markers from marker triads with acceptable error
%ERCOR_fillgaps_RidgidBody(s,'EMGM','RIGIDBODY_PARTIAL_RECONSTRUCTION');


% PLOT simple check for rigidbody errors
figure(hfig); 
sp(4) = subplot2(6,4,5:6,1:4); 
pSE(s,hfig);  % PLOT simple check for rigidbody errors
linkaxes(sp(3:4),'xy');



% CREATE MTATrial object for full session
QuickTrialSetup(S(i),'overwrite',true); 
%QuickTrialSetup(s,'overwrite',true); 
Trial = MTATrial.validate(s.filebase);



if trb,
% ESTIMATE anatomical head position    
    transform_rigidBody(s,false,true);    
    xyz = s.load('xyz','trb');
else
    xyz = s.load('xyz');
    xyz.label = 'trb';
    xyz.key  = 't';
    xyz.name = 'tr_corrected_head';
    xyz.updateFilename(s);    
    xyz.save;
end

% COMPUTE spine interpolated spine
ss = fet_spline_spine(s,'3dssh',xyz,[],true);



% PLOT Session creation report
if report
    hfig = figure(gen_figure_id),
    sp(1) = subplot2(6,4,1:2,1:2);  pZ(Trial);        % PLOT timeseries for z axis
    sp(2) = subplot2(6,4,1:2,3:4);  pXY(Trial);       % PLOT position within xy plane
    sp(3) = subplot2(6,4,3:4,1:4);  pSE(Trial,hfig);  % PLOT simple check for rigidbody errors
end

% LABEL bhv and auxilery behaviors
labelBhv_NN(Trial,'NN0317','jg05-20120317.cof.all','hand_labeled_rev3_jg',...
                  'states',{'walk','rear','turn','pause','groom','sit'});
stc = Trial.load('stc','NN0317');
label_bhv_aux(Trial,stc,[],[],true);

% LABEL reduced bhv and auxilery behaviors
labelBhv_NN(Trial,'NN0317R','jg05-20120317.cof.all','hl_3_jg_r',...
                  'states',{'loc','rear','pause','groom','sit'});
stc = Trial.load('stc','NN0317R');
label_bhv_aux_reduced(Trial,stc,[],[],true);



if nlx,
% COMPUTE neuron quality
    NeuronQuality(s,[],[],[],true); 
% LABEL Theta periods    
    thetaChan = 16;
    Trial = labelTheta(Trial,[],thetaChan,true);
% COMPUTE 3d place fields for the theta state
    pfs_3d_theta(Trial,1,1);
% COMPUTE 2d place field for each state
    pfs_2d_states(Trial);
% COMPUTE 2d state place fields with half sampling bootstrap
    batch_compute_pfstats_bs(Trial);    
    batch_compute_pfstats_bs(Trial,'NN0317',...
                             {'walk','rear','turn','pause',...
                              'lwalk','hwalk','lpause','hpause'});
end
