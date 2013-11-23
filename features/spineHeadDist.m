function feature = spineHeadDist(Session,method,varargin)
[refPair,markerSet] = DefaultArgs(varargin,{Session.Model.gmi({'spine_lower','pelvis_root'}),Session.Model.gmi({'spine_middle','spine_upper','head_back','head_front'})});

Session = Session.load_ang;

feature = Session.ang(:,4,5,3);
