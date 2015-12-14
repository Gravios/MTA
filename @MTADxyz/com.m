function center_of_mass = com(Data,Model,varargin)
%center_of_mass = com(Session,Model)
%
%Model - MTAModel: MTA object holding marker information
%
%Examples:
%  Find the center of mass of the session model
%    center_of_mass = Session.com(Session.xyz.model);
%
%  Select a model based on a subset of markers from a larger model
%    Model = Session.model.rb({'head_back','head_left','head_front','head_right'});
%    center_of_mass = Session.com(Model);
%   
%
[dim] = DefaultArgs(varargin,{[1:3]});
center_of_mass = mean(Data.subsref(substruct('()',{':',Model.ml(),dim})),2);
end
