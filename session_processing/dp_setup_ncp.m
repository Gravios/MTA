%% MTASession Setup - er02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overwrite = false;
nlx = true;
S = get_session_list('ncp');

i=1;
%thetaChan = [71,71,69,72,73,74,71];


if overwrite, QuickSessionSetup(S(i)); end
s = MTASession.validate(S(i)); 

% IF you forgot to include a vsk upon construction
s.load('xyz');
s.xyz.model =  MTAModel(fullfile(s.spath, [s.name '-' s.maze.name '.vsk']),'-vsk');
s.xyz.save

% PLOT Session creation report
hfig = figure(gen_figure_id),
sp(1) = subplot2(6,4,1:2,1:2);  pZ(s);drawnow;        % PLOT timeseries for z axis
sp(2) = subplot2(6,4,1:2,3:4);  pXY(s);drawnow;       % PLOT position within xy plane
sp(3) = subplot2(6,4,3:4,1:4);  pSE(s,hfig);drawnow;  % PLOT simple check for rigidbody errors

% CORRECT marker swaps 
headMarkers = regexpi(s.xyz.model.ml,'^head_[A-Za-z]*','match');
headMarkers = cellfun(@(x) x, headMarkers(~cellfun('isempty',headMarkers)));
% automatic
%ERCOR_fillgaps_RidgidBody(s,'EMGM','BEST_SWAP_PERMUTATION',[],headMarkers);
% manual
ERCOR_fillgaps_RidgidBody(s,'MANUAL','BEST_SWAP_PERMUTATION',[],headMarkers);

% RECONSTRUCT markers from marker triads with acceptable error
%ERCOR_fillgaps_RidgidBody(s,'EMGM','RIGIDBODY_PARTIAL_RECONSTRUCTION');


% PLOT simple check for rigidbody errors
figure(hfig); 
sp(4) = subplot2(6,4,5:6,1:4);  pSE(s,hfig);  % PLOT simple check for rigidbody errors
linkaxes(sp(3:4),'xy');


% ESTIMATE anatomical head position
transform_rigidBody(s,false,true);

% COMPUTE anatomical head position
xyz = s.load('xyz');
ss = fet_spline_spine(s,'3dssh',xyz,[],true);


%% Create All Trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CREATE MTATrial object for full session
QuickTrialSetup(s,'overwrite',true); 
Trial = MTATrial.validate(s.filebase);


% PLOT Session creation report
hfig = figure(gen_figure_id),
sp(1) = subplot2(6,4,1:2,1:2);  pZ(Trial);        % PLOT timeseries for z axis
sp(2) = subplot2(6,4,1:2,3:4);  pXY(Trial);       % PLOT position within xy plane
sp(3) = subplot2(6,4,3:4,1:4);  pSE(Trial,hfig);  % PLOT simple check for rigidbody errors


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
    Trial = labelTheta(Trial,[],thetaChan(i),true);
% COMPUTE 3d place fields for the theta state
    pfs_3d_theta(Trial,1,1);
% COMPUTE 2d state place fields with half sampling bootstrap
    batch_compute_pfstats_bs(Trial);    
    batch_compute_pfstats_bs(Trial,'NN0317',...
                             {'walk','rear','turn','pause',...
                              'lwalk','hwalk','lpause','hpause'});
end

