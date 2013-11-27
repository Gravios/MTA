function [head_state,head_feature] = head(Session,method,varargin)

switch method

  case 'vel'
    [headThresh,minimum_interval,display] = DefaultArgs(varargin,{2,40,0});
    % Head - speed threshold
    % Feature - spine_lower 
    fet = xyzvel(Session.xyz,gmi(Session.Model.markers,'head_front'),[1 2])*Session.xyzSampleRate/10;
    head_feature = medfilt1(fet,20);
    head = zeros(size(head_feature));
    head(head_feature>headThresh) = 1;
    headPeriods = ThreshCross(head,0.5,minimum_interval);

  case 'com'
    [headThresh,minimum_interval,display] = DefaultArgs(varargin,{2,40,0});
    % Head - speed threshold
    % Feature - center of mass of the head markers
    b4 = Session.Model.rb({'head_back','head_left','head_front','head_right'});
    fet = Session.com(b4);
    fetv = sqrt(sum(diff(fet(:,:,[1 2]),1).^2,3));
    head_feature = medfilt1(fetv,32)*Session.xyzSampleRate/10;
    head = zeros(size(head_feature));
    head(head_feature>headThresh) =1;
    headPeriods = ThreshCross(head,0.5,minimum_interval);

  otherwise
    error(['Method: ' method ' does not exist'])
end

head_state = headPeriods;