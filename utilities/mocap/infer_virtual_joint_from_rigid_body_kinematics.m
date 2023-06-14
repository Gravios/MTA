function infer_virtual_joint_from_rigid_body_kinematics(Session,varargin);
%function transform_rigidBody(Session,varargin);
%
% Estimates the position of the neck by locating the position relative to the head wihch
% minimizes translations of head movements.
%
% Perhaps overly complicated, and inefficient, but it works...
%
% varargin:
%   display       (logical)               false,    display markers and gradients
%   overwrite     (logical)               false,    overwrite gradients
%   scoreFunction (struct)                struct('fun',@mad,'args',{{1}}), 
%   periods       (numeric Nx2)           [],      restrict periods to active movements 
%


% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('display',               false,                                                 ...
                 'overwrite',             false,                                                 ...
                 'scoreFunction',         struct('fun',@mad,'args',{{1}}),                       ... terrible name
                 'periods',               [],                                                    ...
                 'rigidBodyMarkers',      {{'head_back','head_left','head_front','head_right'}}  ...
);
[display, overwrite, scoreFunction, periods, rigidBodyMarkers] =                                 ...
    DefaultArgs( varargin, defargs, '--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

%% LOAD session 
if ischar(Session),
% LOAD session
    Session = MTASession.validate(Session);
else
% !!! CHANGE THIS STUPID BEHAVIOR
% CONVERT trial to session
    trialPeriodsByDefault = false;
    if isempty(periods), 
        periods = Session.sync.copy();         
        trialPeriodByDefault = true;
    end    
    Session = MTASession.validate([Session.name,'.',Session.maze.name,'.all']);
    if trialPeriodsByDefault,         
        periods.data = periods.data-Session.xyz.origin;
        periods.resample(Session.xyz.sampleRate);
        periods = periods.data;
    end
end

% LOAD xyz
% ADD lowpass filtered marker of the head's center of mass
xyz = resample(Session.load('xyz'),20);                       % MTADxyz 

[rigidBodyBasis, headCenterOfMass] = compute_rigidbody_basis_timeseries(xyz,rigidBodyMarkers);
headCenterOfMass = sq(headCenterOfMass);

hcomSpeed = sqrt(sum(clip(multiprod(rigidBodyBasis,circshift(headCenterOfMass,-1)-circshift(headCenterOfMass,1),[2,3],[2]),-100,100).^2,2,'omitnan'));


% GET "appropriate" subset of of data for fitting
nind = (hcomSpeed > 3)  &  (hcomSpeed < 100)  &  (headCenterOfMass(:,3) < 150);


virtualJointCenter = [];
for dim = 1:size(xyz,3)
    offset = [0,0,0];
    offset(dim) = 1;
    grad = [];         
    for shift = [-100:100],
        scom = headCenterOfMass + multiprod(rigidBodyBasis,offset*shift,[2,3],2);
        dcom = circshift(scom,-1)-circshift(scom,1);
        scomVelocity = sum(clip(multiprod(rigidBodyBasis(nind,:,:),dcom(nind,:),[2,3],[2]),-100,100).^2,'omitnan');
        grad = cat(1,grad,log10(scomVelocity));
    end
    [~,mind] = min(grad);
    virtualJointCenter(dim) = mean(zShifts(mind));
end


% LOAD original data
xyz = Session.load('xyz');
[rigidBodyBasis, headCenterOfMass] = compute_rigidbody_basis_timeseries(xyz,rigidBodyMarkers);

%multiprod(rigidBodyBasis,offset*shift,[2,3],2);



figure();
hold('on');
plot(zShifts,grad(:,1),'r')
plot(zShifts,grad(:,2),'b')
plot(zShifts,grad(:,3),'g')

