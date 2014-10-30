function Trial = labelTheta(Trial,varargin)
[Stc,thetaChan,overwrite] = DefaultArgs(varargin,{[],1,false});

if isempty(Stc),
    Stc = Trial.stc.copy;
end


sempty = isempty(Stc{'t'});
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('t')) = [];
    end
    cpath = pwd;
    cd(Trial.spath);
    CheckEegStates(Trial.name,[],[],[],thetaChan,[],'compute',overwrite);
    CheckEegStates(Trial.name,[],[],[],thetaChan,[],'display',false);
    Stc.states{end+1} = theta(Trial);
end

Stc.save(1);
Trial.stc = Stc;
Trial.save;