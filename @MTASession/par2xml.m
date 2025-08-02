function parameters = par2xml(Sesssion,varargin)

par.hash = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz';

if strcmp(par.hash,DataHash(rmfield(par.paramteres))),
    warning('MTA:MTASession:par2xml:ExitOnNoChangesDetected');
    return;
else
% T
    Sesssion.load('nq');
    unitsType = false([size(Sesssion.spk.map,1),1]);
    unitsType(select_units(Sesssion,'pyr','all')) = true;
    uTypes = {'int','pyr'};


    for f = fieldnames(par.parameters.generalInfo)',
        par.parameters.generalInfo.(f{1}).Text = par.parameters.generalInfo.(f{1});
    end



    for f = fieldnames(par.parameters.acquisitionSystem)',
        par.parameters.acquisitionSystem.(f{1}).Text = num2str(par.parameters.acquisitionSystem.(f{1}));
    end



    par.parameters.fieldPotentials.lfpSamplingRate.Text = num2str(par.parameters.fieldPotentials.lfpSamplingRate);


    % IGNORE par.parameters.files
    % IGNORE par.parameters.anatomicalDescription
    % IGNORE par.parameters.spikeDetection

    %--- correct to here ---%
    cind = []

    for u = 1:numel(par.parameters.units)
        units.unit{u}.group.Text       = num2str(par.parameters.units(u).group);
        units.unit{u}.cluster.Text     = num2str(par.parameters.units(u).cluster);
        units.unit{u}.structure.Text   = par.parameters.units(u).structure;
        units.unit{u}.type.Text        = par.parameters.units(u).type;    
        units.unit{u}.isolationDistance.Text = num2str(par.parameters.units(u).isolationDistance);
        units.unit{u}.quality.Text = num2str(par.parameters.units(u).quality);
        units.unit{u}.notes.Text        = par.parameters.units(u).notes;    
    end    

    par.parameters.units = units;

    par.hash = DataHash(par.paramteres);
end