function [rigidBodyBasis, hcom] = compute_rigidbody_basis_timeseries(xyz,varargin)
% function [rigidBodyBasis, hcom] = compute_rigidbody_basis_timeseries(xyz,varargin)
%
% computes the basis of the rigidbody specified by the markers centered on the mean 
% position of the markers for all timepoints in the xyz series.
%
% varargin:
%   markers (cellstr 1xN)   {'head_back','head_left','head_front','head_right'}, names of markers
%

% DEFARGS ------------------------------------------------------------------------------------------
defargs   = struct('markers', {'head_back','head_left','head_front','head_right'});

[markers] = DefaultArgs( varargin, defargs, '--struct');
%---------------------------------------------------------------------------------------------------

ERR_PREFIX = 'MTA:utilities:mocap:compute_rigidbody_basis_timeseries:';

assert( isa(xyz,'MTADxyz'),           [ ERR_PREFIX, ' xyz is not an MTADxyz object']);
assert( numel(markers)>1  ,           [ ERR_PREFIX, ' the number of markers must be greater than 1']);
assert( all(cellfun(@ischar,markers)),[ ERR_PREFIX, ' one or marker names are non-char type']);

% GENERATE MTAModel object for the set of markers within the target rigid body
rigidBodyModel = xyz.model.rb(markers);     % MTAModel

% COMPUTE the center of mass of the rigid body ∀ timepoints
hcom = xyz.com(rigidBodyModel);

% GENERATE orthogonal basis, origin: head's center of mass
nz = -cross(xyz(:,markers{1},:)-hcom,xyz(:,markers{2},:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 

% GENERATE orthogonal basis, origin: head's center of mass
ny = cross(nz,xyz(:,markers{1},:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));

% GENERATE orthogonal basis, origin: head's center of mass
nx = cross(nz,ny);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));

rigidBodyBasis = cat(3,sq(nx),sq(ny),sq(nz));



