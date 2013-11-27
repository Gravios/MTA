function crtscpt(Session,varargin)
[modelIndex,display,markerSubset,depth] = DefaultArgs(varargin,{[],1,{},5});

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end

if isempty(markerSubset),
    for i = 1:length(Session.Model.ml()),
        if regexp(Session.Model.Markers{i}.name,'^head'),
            markerSubset{end+1} = Session.Model.Markers{i}.name;
        end
    end
end

if display,
    figure
    hbflr = Session.transformOrigin('head_back','head_front',{'head_left','head_right'});
    hrlbf = Session.transformOrigin('head_right','head_left',{'head_back','head_front'});
    sp1 = subplot(211);
    plot([hbflr.transVec(:,:,2),hrlbf.transVec(:,:,2)])
end

Session = Session.correctRigidBody(markerSubset,depth,modelIndex);
Session = Session.correctPointErrors(markerSubset);
Session.save();

if display,
    hbflr = Session.transformOrigin('head_back','head_front',{'head_left','head_right'});
    hrlbf = Session.transformOrigin('head_right','head_left',{'head_back','head_front'});
    sp2 = subplot(212);
    plot([hbflr.transVec(:,:,2),hrlbf.transVec(:,:,2)])
    linkaxes([sp1,sp2],'xy')
end

