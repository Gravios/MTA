
%% 2015 07 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20150728
%    create video of rear 



%% 2015 08 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% req20150914
%    Run tsne on a session and overlay the hand labels, and then plot
%    mean z-scores for each feature over the tsne space. 

% req20150915
%    precursor to the final version which is req20150914.m

% req20150925
%    Run tsne on a session and overlay the hand labels, and then plot
%    mean z-scores for each feature over the tsne space. 



%% 2015 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20151007
%    logistic regression models 
%    label other sessions with hand labeled states, using the lgr model
%    trained on the data set from jg05-20120317, and comapare with
%    hand labels.

% req20151009
%    All JPDFs for fet_tsne over all sessions


%% 2015 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% req20160114
%    tau masked boundaries for confusion matrix computation

% req20160126
%    gausshmm segmentation of features into states

% req20160127
%    head vs body vel JPDF w/ state contours

% req20160128
%    Train nn classifiers with bhv_nn_multi_patternnet


%% 2016 02 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20160212
%    playing with feature copulas 

% req20160224
%    precursor to req20160310


%% 2016 03 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20160310
%    miha feature selection

% req20160315
%    preliminary analysis of optogenetic rat data and mouse data
%    with optitrack


%% 2016 04 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% req20160411
%    Examine inter-session overlap in t-sne space

