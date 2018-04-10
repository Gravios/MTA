
%% 2015 07 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20150728
%    create video of rear 



%% 2015 08 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20150817
%    Feature exploration and visualization 

% req20150821 
%    status: inactive
%    REQ 1.0 Video Example of skeleton and feature phase space
%    REQ 2.0 Err,... some more fanangling with walking features
%    REQ 3.1 Rear Body pitch vs diff(USpitch)
%    REQ 3.2 Rear US pitch vs US log10(abs(diff(ang)))
%    REQ 4.1 Walking Segmentation



%% 2015 09 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20150914
%    Run tsne on a session and overlay the hand labels, and then plot
%    mean z-scores for each feature over the tsne space. 

% req20150915
%    precursor to the final version which is req20150914.m

% req20150925
%    Run tsne on a session and overlay the hand labels, and then plot
%    mean z-scores for each feature over the tsne space. 



%% 2015 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20151007
%    logistic regression models 
%    label other sessions with hand labeled states, using the lgr model
%    trained on the data set from jg05-20120317, and comapare with
%    hand labels.

% req20151009
%    All JPDFs for fet_tsne over all sessions


%% 2015 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20151103
%    Neural network stuff
%    1. Train a neural network on the data of jg05-20120317. 
%    2. Behavior of other animals is then classified with the network. 
%    3. The Results are compared to the hand labeled data with a confussion matrix

% req20151118.m
%    Examination of rats' casting motions during exploration
%    1. Identify types of head casting
%    2. Characterize phiysical characteristics of each casting type
%    3. Determine relationship between casting and sniffing

% req20151123.m
%    Primary Goals
%    Feature correction and normalization between subjects



%% 2015 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20151203:
%    tSNE over all jg05 sessions excluding jg05-20120317
%    jg05-20120317 features are then mapped onto tsne space
%
%        alt1:
%        Equal bootstrap noisy trim
%        


% req20151204
%    tsne mapping of spike on to 2D 

% req20151216
%    neural network training validation randomly selected data
%    points.
%
%    TODO: random behavior epoch sampling.
%


% req20151217
% General testing of the bhv_nn wraper for the neural network
% packages of matlab


%% 2016 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20160114
%    tau masked boundaries for confusion matrix computation

% req20160126
%    gausshmm segmentation of features into states

% req20160127
%    head vs body vel JPDF w/ state contours

% req20160128
%    Train nn classifiers with bhv_nn_multi_patternnet


%% 2016 02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20160212
%    playing with feature copulas 

% req20160224
%    precursor to req20160310


%% 2016 03 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20160310
%    miha feature selection

% req20160315 (stasis)
%
%  Final_Forms:
%
%  Description:
%    preliminary analysis of optogenetic rat data and mouse data
%    with optitrack


%% 2016 04 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20160411
%    Examine inter-session overlap in t-sne space

% req20160419
%    Testing Fast Artificial Neural Network (FANN) c library



%% 2016 05 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20160510
%    Report feature distributions for a set of mappings and
%    normalizations

% req20160517
%    find plane through lower 4 head markers

% req20160518
%    Accumulate stats for all label/train pairings for
%    bhv_multi_patternnet




%% 2016 06 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.1

% req20160615 (stasis) ---------------------------------------------
%    Check the transform_rigidbody script against emperical values

% req20160620 (stasis) ---------------------------------------------
%    Examination of head pitch dependence on joint estimation from
%    rigidbody dynamics

% req20160621 (statis) ---------------------------------------------
%    Examination of translated vs raw data rhm coherence with ncp





%% 2016 07 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.2

% req20160719 ----------------------------------------------------
%
%  Status: retired
%  Type: Analysis
%  Final_Forms: 
%    analysis/exp_teleport.m
%    analysis/exp_teleport_0823.m
%  Description:
%    Initial script for the shift teleportation experiment of the
%    VR system. 
%  Bugs: NA



% req20160721 ----------------------------------------------------
%
%  Status: Active
%  Type: Analysis
%  Final_Forms: NA
%  Description: 
%    Exploration of head features for the segmentation of states
%    which have previously required body features (postures/kinematics)
%  Bugs: NA



% req20160730 ----------------------------------------------------
%
%  Status: retired
%  Type: Utility
%  Final_Forms: utilities/orient_head_to_rhm.m
%  Description:
%    Create script which from a rigid body, will find a vector in
%    the direction of locomotion which best represents the
%    orientation of the head. Uses rhythmic head motion as primary
%    classification feature. Assumes the front to back to front
%    motion of the head during respiration is relatively contstant
%    between animals.
%  Bugs: only works for Nick's data structure



%% 2016 08 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.2


% req20160815 ----------------------------------------------------
%
%  Title: Organization of rhythmic head motion and gait 
%  Type: Analysis
%  Status: retired
%  Final_Forms: NA
%  Description:
%    Script for subjective evaluation of MTAAknnpfs_bs parameters
%  Bugs: NA

% req20160822 ----------------------------------------------------
%
%  Title: Data processing for session ER06-20130612
%  Type: Utility
%  Status: retired
%  Final_Forms: NA
%  Description:
%    Abandoned code for processing a session.
%  Bugs: NA

% req20160831 ----------------------------------------------------
%
%  Title: Train neural networks with alternate marker models
%  Type: Utility
%  Status: active
%  Final_Forms: NA
%  Description:
%    Train neural networks with non-standard marker models. B4H5
%    has been the standard model for all previous work.
%  Bugs:
%    Requires a feature script which generates all features
%    available for the given model




%% 2016 09 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.2


% req20160901 ----------------------------------------------------
%
%  Status: retired
%  Type: Analysis
%  Final_Forms: NA
%  Description:
%    Unknow purpose
%  Bugs: NA

% req20160927 ----------------------------------------------------
%
%  Status: retired
%  Type: Analysis
%  Final_Forms: NA
%  Description:
%    Explore use of k-means for feature clustering
%  Bugs: NA

% req20160929 ----------------------------------------------------
%
%  Status: retired
%  Type: Utility
%  Final_Forms: NA
%  Description:
%    Script for subjective evaluation of MTAAknnpfs_bs parameters
%  Bugs: NA



%% 2016 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.2


% req20161006 ----------------------------------------------------
%
%  Status: retired
%  Type: Analysis 
%  Project: MjgEd2016
%  Final_Forms: NA
%  Description:
%    Exploring copulas
%  Bugs: NA



%% 2016 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.2


% req20161103 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Project: MjgER2016
%  Final_Forms: NA
%  Description:
%    checking MTAAknnpfs shuffle feature for thresholding during
%    place field statistics calculated from threshholded patches
%  Bugs: NA

% req20161116 ----------------------------------------------------
%
%  Status: active
%  Type: Utility
%  Project: MTA
%  Final_Forms: NA
%  Description:
%    optimize points in head rigidbody based on other markers
%  Bugs: NA

% req20161124 ----------------------------------------------------
%
%  Status: active
%  Type: Utility
%  Project: MTA
%  Final_Forms: NA
%  Description:
%    some sort of spike sorting probe
%  Bugs: NA

% req20161125 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis
%  Project: MjgEd2016
%  Final_Forms: NA
%  Description:
%    playing with some spectral features 
%  Bugs: NA



%% 2016 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% format: v0.2

% req20161206 ----------------------------------------------------
%
%  Status: active
%  Type: analysis
%  Final_Forms: NA
%  Description:
%    Examining data sets with hip markers, and walking features.
%  Bugs: NA

% req20161212 ----------------------------------------------------
%
%  Status: retired
%  Type: Utility 
%  Final_Forms: convert_stc2sts.m
%  Description:
%    convert Stc states into 1250Hz sts files for general lab 
%    analysis
%  Bugs: NA

% req20161230 ----------------------------------------------------
%
%  Status: active
%  Type: Utility 
%  Final_Forms: NA
%  Description:
%    Experimentation with visualization of 3d place fields. Using
%    surf to display xy slices along the z-axis.
%  Bugs: NA

% req20161231 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description:
%    state transition unit ccgs
%  Bugs: NA




%% 2017 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2017 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20170126 ----------------------------------------------------
%
%  Status: active
%  Project: MjgEd2016
%  Type: Utility 
%  Final_Forms: NA
%  Description:
%    Examine current inter-session/subject feature coregistration
%  Bugs: NA


% req20170127 ----------------------------------------------------
%
%  Status: active
%  Project: MjgEd2016
%  Type: Utility 
%  Final_Forms: NA
%  Description:
%    Segment locomotive trajectories into in/out bound from 
%    homebase locations
%  Bugs: NA


%% 2017 02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% req20170201 ----------------------------------------------------
%
%  Status: active
%  Project: Mjgvk2017
%  Type: Utility 
%  Final_Forms: NA
%  Description: Convert MTA objects to 2d mat for vassia colaboration
%  Bugs: NA

% req20170211 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Check out neural network reconstruction of timeseries
%  Bugs: NA

% req20170217 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Labeling statistics 
%               Req 1.0: stats between neural network classifier (NNC) labels
%               Req 1.1: stats between expert labeler classifier (ELC) labels 
%               Req 1.2: stats between expert labeler classifier (ELC) labels 
%                        and neural network classifier (NNC) labels 
%  Bugs: NA
 



%% 2017 06 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% req20170605
%
%  Status: inactive
%  Type: Analysis 
%  Final_Forms: fet_rhmPCA
%  Description: comparison 
%  Bugs: NA
%



%% 2017 07 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20170707 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Behavior subtyping with t-sne 
%  Bugs: NA

% req20170715 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: load_normalization_mapminmax.m
%  Description: multi session normalization
%  Bugs: NA

% req20170726 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Quick ccg(behavior,unit) in the cortex
%  Bugs: NA
 
% req20170728 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Labeling statistics 
%  Bugs: NA
 
% req20170731 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Labeling statistics 
%  Bugs: NA




%% 2017 08 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20170802 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Behavior subtyping with t-sne 
%  Bugs: NA


% req20170811 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Behavior subtyping with t-sne 
%  Bugs: NA


%% 2017 09 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20170901 ----------------------------------------------------
%  Status: inactive
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Find head only features for grooming detection
%  Bugs: NA

% req20170906 ----------------------------------------------------
%  Status: inactive
%  Type: Analysis 
%  Final_Forms: NA
%  Description: exploration of other rhythmic head motions
%  Bugs: NA

% req20170912 ----------------------------------------------------
%  Type: Analysis (distraction)
%  Final_Forms: NA
%  Description: behavior labeling with recursive neural network 
%  Bugs: NA

% req20170913 ----------------------------------------------------
%  Type: Analysis
%  Final_Forms: NA
%  Description: Single Step Analysis
%  Bugs: NA

% req20170914 ----------------------------------------------------
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: revisiting repair rigidbody partial reconstruction 
%  Bugs: NA

% req20170921 ----------------------------------------------------
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Head body pitch azimuth xcorr
%  Bugs: NA

% req20170924 ----------------------------------------------------
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Hippocampal Layer depth based on ephys data
%  Bugs: NA

% req20170925 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: new rhm feature base on proper filtering
%  Bugs: NA

% req20170926 ----------------------------------------------------
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: examination of new RectFilter and usefulness in fet_rhm
%  Bugs: NA

%% 2017 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20171020 ----------------------------------------------------
%  Status: retired
%  Type: Analysis
%  Final_Forms: NA
%  Description: exploration of other rhythmic head motions
%  Bugs: NA

% req20171023 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Second quick look into lfp stuff
%  Bugs: NA

% req20171026 ----------------------------------------------------
%  Status: active
%  Type: Test
%  Final_Forms: NA
%  Description: Debug and Characterize the label_bhv_reduced_aux classifier
%  Bugs: NA

% req20171031 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: MjgER2016_figure3.m
%  Description: processing and display of placefields 
%  Bugs: NA




%% 2017 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20171101 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: MjgER2016_drzfields.m
%  Description: directional rate zone (DRZ) versus various features
%  Bugs: NA


% req20171107 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: State occupancy breakdown
%  Bugs: NA



%% 2017 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20171202 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Optimizing place field and drz parameters
%  Bugs: NA

% req20171205 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Determine the respiration frequency distribution during 
%               immobility
%  Bugs: NA

% req20171212 ----------------------------------------------------
%  Status: active
%  Type: ProcessData
%  Final_Forms: NA
%  Description: Go through Eduardo's rats which contain pressure 
%               sensors to collect statistics.
%  Bugs: NA

% req20171213 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: drzfield population stats
%  Bugs: NA

% req20171214 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Place fields in non theta periods minus ripples
%  Bugs: NA

% req20171215 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: use placefield rate and area to search for point 
%               spatial representation
%  Bugs: NA

% req20171216 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: parameter scripts for testing locate_rhm_ncp_coherence_point
%  Bugs: NA

% req20171220 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Phase precession analysis
%  Bugs: NA




%% 2018 01 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20180109 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Refine place field selection, based on spatial information
%               of each state
%  Bugs: NA

% req20180117 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Quick Gamma Burst Analysis
%  Bugs: NA

% req20180120 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: testing new formulation of rhm
%  Bugs: NA

% req20180123 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: RateMap(Pitch,Height) where drz := [-0.5,0.5]
%  Bugs: NA

% req20180127 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Comparison of trb hcom and nose prediction for placefields
%  Bugs: NA

% req20180129 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: M: spatial self-representation
%  Bugs: NA



%% 2018 02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20180206 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: check feature coregistration quality of xyz preprocessing
%               algorithm modification SPLINE_SPINE_HEAD_EQD option
%  Bugs: NA

% req20180214 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: PlotPF.m:var:OccThresh
%  Keywords: Placefields, PlotPF, Threshold
%  Description: establish new criteria for placefield occupacy 
%               threshold
%  Bugs: Temporary threshold established.
%        Replace with mean velocity map corrections.

% req20180216 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Keywords: Placefields, drz, behavioral transition
%  Description: Computing the mean spike theta phase as a function
%               of drz and behavior onset/offset
%  Bugs: NA

% req20180222 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Keywords: Placefields, drz, mds, t-sne
%  Description: Compute rate map of reduced dimesion subspace
%               see fet_rds.m
%  Bugs: NA


% req20180227 ----------------------------------------------------
%  Status: active
%  Type: Diagnostic
%  Final_Forms: NA
%  Keywords: Placefields, sampleRate, bootstrap
%  Description: placefield estimation as function of xyz samplerate
%  Bugs: NA
%  Notes: Reduced sample rate (16Hz) works sufficiently well and 
%         saves computational resources 




%% 2018 03 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20180301 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: lfp analysis cycle by cycle comparison
%               
%  Bugs: NA

% req20180309 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Distance vs drz restricted behavior field egenvector scores
%               
%  Bugs: NA

% req20180319 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Compute z-score of randomly shuffled spike time chunks 
%               projected on each eigen vector, as computed by 
%               req20180123_ver5.m
%  Bugs: NA

% req20180322 ----------------------------------------------------
%  Status: active
%  Type: Utility
%  Final_Forms: NA
%  Description: Review unit placefields 
%  Bugs: NA

% req20180323 ----------------------------------------------------
%  Status: active
%  Type: Utility
%  Final_Forms: NA
%  Description: inter-session unit corregistration
%  Bugs: NA




% req20180123_ver5.m 
% compute and decompose unit rate multidimensional tunning

% req20180319.m
% Compute z-score of randomly shuffled spike time chunks 
% projected on each eigen vector, as computed by 
% req20180123_ver5.m

% req20180322.m
% Review unit placefields

% req20180323.m
% inter-session unit corregistration