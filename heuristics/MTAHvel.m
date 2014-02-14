function [data,sampleRate,label,key] = MTAHvel(Session,varargin)

    [threshold,label,key] = DefaultArgs(varargin,{2,'vel','v'});
    
    txyz = Session.xyz.copy;
    
    Session.xyz.clear;
    Session.xyz.load(Session);

    order = round(.25*Session.xyz.sampleRate);
    if mod(order,2)==0,order = order+1;end

    Session.xyz.filter(gausswin(order)./sum(gausswin(order)));

    data = ThreshCross(Session.vel(Session.trackingMarker,[1,2]),threshold,10);
    sync = Session.xyz.sync.copy;
    
    data = MTADepoch(Session.spath,...
                     Session.filebase,...
                     data,...
                     Session.xyz.sampleRate,...
                     sync,...
                     sync(1),...
                     label,key);

    Session.xyz = txyz;
end
