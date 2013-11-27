function feature = spineAngZ(Session,method,varargin)
[refPair,markerSet] = DefaultArgs(varargin,{Session.Model.gmi({'spine_lower','pelvis_root'}),Session.Model.gmi({'spine_middle','spine_upper','head_back','head_front'})});

Session = Session.load_ang;

feature = unwrap(sq(Session.ang(:,refPair(1),markerSet,2))-repmat(sq(Session.ang(:,refPair(1),refPair(2),2)),1,length(markerSet)));
