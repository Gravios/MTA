function Data = subset(Data,markers)
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
Data.data = Data.data(:,Data.model.gmi(markers),:);
Data.model = Data.model.rb(markers);
end
