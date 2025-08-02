function parameters = xml2parameters(Trial)

[par] = xml2struct(fullfile(Trial.spath,Trial.name));
parameters = par.parameters;
par.hash = 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz';

Trial.load('nq');
unitsType = false([size(Trial.spk.map,1),1]);
unitsType(select_units(Trial,'pyr','all')) = true;
uTypes = {'int','pyr'};


for f = fieldnames(par.parameters.generalInfo)',
    par.parameters.generalInfo.(f{1}) = par.parameters.generalInfo.(f{1}).Text;
end

for f = fieldnames(par.parameters.acquisitionSystem)',
    par.parameters.acquisitionSystem.(f{1}) = str2num(par.parameters.acquisitionSystem.(f{1}).Text);
end

par.parameters.fieldPotentials.lfpSamplingRate = str2num(par.parameters.fieldPotentials.lfpSamplingRate.Text);


% IGNORE par.parameters.files
% IGNORE par.parameters.anatomicalDescription
% IGNORE par.parameters.spikeDetection

cind = []

for u = 1:numel(par.parameters.units.unit)

    groupId = str2num(par.parameters.units.unit{u}.group.Text);
    clusterId = str2num(par.parameters.units.unit{u}.cluster.Text);

    cind(u) = find(ismember(Trial.spk.map(:,[2,3]),[groupId,clusterId],'rows'));
    
    units(u) = struct('group',                    groupId,                                    ...
                      'cluster',                  clusterId,                                  ...
                      'structure',                par.parameters.units.unit{u}.structure.Text,...
                      'type',                     uTypes{unitsType(cind(u))+1},               ...
                      'isolationDistance',        Trial.nq.eDist(cind(u)),                    ...
                      'quality',                  1,                                          ...
                      'notes',                    par.parameters.units.unit{1}.notes.Text     ...
    );
end    
[~,scind] = sort(cind);

par.parameters.units = units(scind);

par.hash = DataHash(rmfield(par.paramteres));
