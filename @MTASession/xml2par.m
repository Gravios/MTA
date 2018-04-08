function Session = xml2par(Session)

[par] = xml2struct(fullfile(Session.spath,[Session.name,'.xml']));
parameters = par.parameters;
clear('par');


Session.load('nq');
unitsType = false([size(Session.spk.map,1),1]);
unitsType(select_units(Session,'pyr','all')) = true;
uTypes = {'int','pyr'};


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

cind = []
for u = 1:numel(parameters.units.unit)

    groupId = str2num(parameters.units.unit{u}.group.Text);
    clusterId = str2num(parameters.units.unit{u}.cluster.Text);

    cind(u) = find(ismember(Session.spk.map(:,[2,3]),[groupId,clusterId],'rows'));
    
    units(u) = struct('group',                    groupId,                                    ...
                      'cluster',                  clusterId,                                  ...
                      'structure',                parameters.units.unit{u}.structure.Text,...
                      'type',                     uTypes{unitsType(cind(u))+1},               ...
                      'isolationDistance',        Session.nq.eDist(cind(u)),                    ...
                      'quality',                  1,                                          ...
                      'notes',                    parameters.units.unit{1}.notes.Text     ...
    );
end    
[~,scind] = sort(cind);
parameters.units = units(scind);


Session.parametersHash = DataHash(parameters);
Session.parameters     = parameters;