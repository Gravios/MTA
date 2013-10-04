function immobile_state = immobile(Session,varargin)
% immobility 
% Feature - not walk / not rear
% Heuristic - exclusion

immobile = ones(size(Session.xyz,1),1);
immobile([0;Session.Bhv.states(]|Session.Bhv.rear) = 0;
immobele(1)=0;
