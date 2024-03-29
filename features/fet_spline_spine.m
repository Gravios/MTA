function ssp = fet_spline_spine(Trial,varargin)
% function ssp = fet_spline_spine(Trial,label,xyzMode,markers,overwrite)
% create an spline interpolated data set for 100 points between tail and head
%
% DEFARGS 
%
%  label:    string,      '3dssh',
%  xyzMode:  string,      'trb',                                               
%  markers:  CellArray,   {{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}},
%  overwrite:Logical,     false                                                
%


% DEFARGS ------------------------------------------------------------------------------------------
%Trial = MTASession.validate(Trial);
Session = MTASession.validate([Trial.name,'.',Trial.maze.name,'.all']);
defargs = struct('label',                   '3dssh',                                             ...
                 'xyzMode',                 'trb',                                               ...
                 'markers',  {{'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'}},...
                 'overwrite',               false                                                ...
);
[label,xyzMode,markers,overwrite] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

xyz = Session.load('xyz',xyzMode);

if isempty(xyz.model.gmi('hcom'))
    xyz.addMarker('hcom',...     Name
                  [.7,0,.7],...  Color
                  {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
                   {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
                   {'head_front','hcom',[0,0,255]},...
                   {'head_right','hcom',[0,0,255]}},... 
                  xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
end

                  
% $$$ switch label,
% $$$   case '3dss'
% $$$     markers = {'spine_lower','pelvis_root','spine_middle','spine_upper'};
% $$$   case '3dssh'
% $$$     markers = {'spine_lower','pelvis_root','spine_middle','spine_upper','hcom'};
% $$$ end

% DELBLOCK if version>=3,
if isempty(Session.fet),
    Session.fet = MTADfet(Trial.spath,...
                        [],...
                        [],...
                        [],...
                        Session.sync.copy,...
                        Session.sync.data(1),...
                        []);                  
end
% DELBLOCK end

filename = cell2mat(cellstr_append_str([Session.spath,'/'],list_files(Session.name,[Session.trialName,'.fet.',label,'.'])));
if overwrite||isempty(filename),
    txyz = xyz(:,markers,:);
    pnts = zeros([xyz.size(1),size(fnplt(cscvn(sq(txyz(find(nniz(xyz),1,'first'),:,:))'))',1),3]);
    for ind = find(nniz(xyz))',
        try
            pnts(ind,:,:) = fnplt(cscvn(sq(txyz(ind,:,:))'))';
        end
    end
    name = '3d spline interpolated spine'; key = 's';
    ssp = MTADfet.encapsulate(Session,...
                              pnts,...
                              xyz.sampleRate,...
                              name,label,key);
    ssp.updateFilename(Session);
    ssp.save;
else
    xyz = Trial.load('xyz',xyzMode);
    ssp = Trial.load('fet',label);
end

ssp.resample(xyz);

% END MAIN -----------------------------------------------------------------------------------------
