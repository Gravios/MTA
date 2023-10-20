function xyz = infer_virtual_joint_from_rigidbody_kinematics(Trial,varargin);
%function transform_rigidBody(Trial,varargin);
%
% Estimates the position of the neck by locating the position relative to the head wihch
% minimizes translations of head movements.
%
% Perhaps overly complicated, and inefficient, but it works...
%
% varargin:
%   sampleRate    (Numeric)               40Hz, sampling rate used during coputations
%   overwrite     (logical)               false,    overwrite xyz object with the crb
%   
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(      'sampleRate',    40,                                                      ... 
                 'rigidBodyMarkers',    {{'head_back','head_left','head_front','head_right'}},   ...
                 'auxiliaryMarkers',    {{'head_top'}},                                          ...
                              'key',    'c',                                                     ...
                            'label',    'crb',                                                   ...
                             'name',    'Corrected Rigid Body',                                  ...
                        'overwrite',    false                                                    ...
);
[sampleRate, rigidBodyMarkers, auxiliaryMarkers, key, label, name, overwrite] = ...
    DefaultArgs( varargin, defargs, '--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

% LOAD session
Session = MTASession.validate(Trial);

xyz = resample(Session.load('xyz'),sampleRate);

[rigidBodyBasis, headCenterOfMass] = compute_rigidbody_basis_timeseries( xyz, rigidBodyMarkers);

headSpeed = sqrt(sum(clip(multiprod(rigidBodyBasis,circshift(sq(headCenterOfMass),-1)-circshift(sq(headCenterOfMass),1),[2,3],[2]).*sampleRate/10/2,-1000,1000).^2,2,'omitnan'));

% GET "appropriate" subset of of data for fitting
nind = ( headSpeed  > 5 )  &  ( headSpeed  < 100 )  &  ( headCenterOfMass(:,3) < 150 ); 
% Units              cm/s                    cm/s                                 mm

shifts = [-100:100]; 
% Units:        mm

% SWITCH from [ N x 1 x 3 ] to [ N x 3 ]
headCenterOfMass = permute(headCenterOfMass,[1,3,2]);

virtualJointCenter = zeros([1,size(xyz,3)]);
for dim = 1:size(xyz,3)
    offset = [0,0,0];
    offset(dim) = 1;
    grad = [];         
    for shift = shifts,
        scom = headCenterOfMass + multiprod(rigidBodyBasis,offset*shift,[2,3],2);
        dcom = circshift(scom,-1)-circshift(scom,1);
        scomVelocity = sum(clip(multiprod(rigidBodyBasis(nind,:,:),dcom(nind,:),[2,3],[2]),-100,100).^2,'omitnan');
        grad = cat(1,grad,log10(scomVelocity));
    end
    [~,mind] = min(grad);
    virtualJointCenter(1,dim) = mean(shifts(mind));
end

% LOAD original data
xyz = Session.load('xyz');
[rigidBodyBasis, headCenterOfMass] = compute_rigidbody_basis_timeseries(xyz,rigidBodyMarkers);

rigidBodyShift = permute(multiprod(rigidBodyBasis,virtualJointCenter,[2,3],2),[1,3,2]);

% ADD shifted headCenterOfMass marker to the MTADxyz object
xyz.addMarker('hcom',                                 ... marker name
              [0.5,1,0.5],                            ... marker color
              { {'head_back', 'hcom',[0,0,1]},        ... stick args
                {'head_front','hcom',[0,0,1]},        ...
                {'head_left', 'hcom',[0,1,0]},        ...
                {'head_right','hcom',[1,0,0]} },      ...
              headCenterOfMass+rigidBodyShift);         % DATA

% SHIFT the rigidbody markers to the virtual joint
for marker = rigidBodyMarkers
    marker = marker{1};
    xyz.data(:,xyz.model.gmi(marker),:) = xyz(:,marker,:) + rigidBodyShift;
end

% SHIFT the auxiliary markers to the virtual joint
auxiliaryMarkers = nonzeros(xyz.model.gmi(auxiliaryMarkers));
if 0 > numel(auxiliaryMarkers),
    for marker = auxiliaryMarkers 
            xyz.data(:,marker,:) = xyz(:,marker,:) + rigidBodyShift;
    end
end

xyz.key  = key;
xyz.label = label;
xyz.name = name;
xyz.updateFilename(Session);      

if overwrite || ~exist(xyz.fpath,'file')
% RENAME & SAVE the MTADxyz object
    xyz.save();
end
