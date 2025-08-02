function Session = update_parameters(Session)
% function Session = update_parameters(Session)
%
% Convert xml parameter file to structure
% 


[par] = xml2struct(fullfile(Session.spath,[Session.name,'.xml']));
parameters = par.parameters;
clear('par');


for f = fieldnames(parameters.generalInfo)',
    parameters.generalInfo.(f{1}) = parameters.generalInfo.(f{1}).Text;
end

for f = fieldnames(parameters.acquisitionSystem)',
    parameters.acquisitionSystem.(f{1}) = str2num(parameters.acquisitionSystem.(f{1}).Text);
end

parameters.fieldPotentials.lfpSamplingRate = str2num(parameters.fieldPotentials.lfpSamplingRate.Text);


% IGNORE parameters.files
% IGNORE parameters.anatomicalDescription
% IGNORE parameters.spikeDetection



Session.load('nq');
unitsType = false([size(Session.spk.map,1),1]);
unitsType(select_units(Session,'pyr','all')) = true;
uTypes = {'int','pyr'};
cind = [];

for u = 1:numel(parameters.units.unit)

    if isfield(parameters.units.unit{u}.group,'Attributes')
        for f = fieldnames(parameters.units(u))',
            f = f{1};
            if strcmpi(parameters.units.unit{u}.(f).Attributes.type,'Numeric')
                units(u).(f) = str2num(parameters.units.unit{u}.(f).Text);
            elseif strcmpi(parameters.units.unit{u}.(f).Attributes.type,'Char')
                units(u).(f) = parameters.units.unit{u}.(f).Text;
            end
        end
    else
% ASSUME xml has never been read by MTA
% CAST each field explicitly 
% REORDER to match MTASpk:map
        units(u).group     = str2num(parameters.units.unit{u}.group.Text);
        units(u).cluster   = str2num(parameters.units.unit{u}.cluster.Text);
        
        cind(u) = find(ismember(Session.spk.map(:,[2,3]),[units(u).group,units(u).cluster],'rows'));
        
        units(u).structure = parameters.units.unit{u}.structure.Text;
        units(u).type      = uTypes{unitsType(cind(u))+1};
        units(u).isolationDistance = round(Session.nq.eDist(cind(u)),4);
        units(u).quality   = 1;
        units(u).notes     = parameters.units.unit{u}.notes.Text;
    end
    
end

if ~isempty(cind),
% REORDER to match MTASpk:map    
    [~,scind] = sort(cind);
    units = units(scind);
end
parameters.units = units;

nq = Session.nq;
nq = rmfield(nq,'AvSpk');
for f = fieldnames(nq)'    
    for u = 1:numel(parameters.units)
        parameters.units(u).(f{1}) = nq.(f{1})(u);
    end
end


Session.parametersHash = DataHash(parameters);
Session.parameters     = parameters;