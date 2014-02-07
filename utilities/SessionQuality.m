function squal = SessionQuality(Session)

Session.load('nq');

squal.fname    = Session.filebase;
squal.nMarkers = Session.model.N;
squal.nClu     = size(Session.spk.map,1);
squal.nCluPyr  = sum(Trial.nq.SpkWidthR>0.4&Trial.nq.TimeSym>2.5);
squal.nCluInt  = sum(~(Trial.nq.SpkWidthR>0.4&Trial.nq.TimeSym>2.5));
squal.nCluG30  = sum(Session.nq.eDist>=30);

%% Total Distance Traveled
dxyz = diff(sq(Session.xyz(:,Session.Model.gmi(Session.trackingMarker),1:2)),2);
xySpeed = sqrt(sum(dxyz.^2,2));
squal.tDistance = sum(xySpeed)/1000;    



%% LIA Time

%% Theta Time

%% Rearing Time

%% Walking Time

%% HWalking Time

%% LWalking Time




