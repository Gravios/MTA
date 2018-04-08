function save_parameters(Session)

if strcmp(Session.parametersHash,DataHash(Session.parameters)),
    warning('MTA:MTASession:par2xml:NoChangesDetected');
    warning('MTA:MTASession:par2xml:AbortSave');    
    return;
else    

    [par] = xml2struct(fullfile(Session.spath,[Session.name,'.xml']));

    for f = fieldnames(Session.parameters.generalInfo)',
        par.parameters.generalInfo.(f{1}).Text = ...
            Session.parameters.generalInfo.(f{1});
    end

    for f = fieldnames(Session.parameters.acquisitionSystem)',
        par.parameters.acquisitionSystem.(f{1}).Text = ...
            num2str(Session.parameters.acquisitionSystem.(f{1}));
    end

    par.parameters.fieldPotentials.lfpSamplingRate.Text = ...
        num2str(Session.parameters.fieldPotentials.lfpSamplingRate);

    % IGNORE Session.parameters.files
    % IGNORE Session.parameters.anatomicalDescription
    % IGNORE Session.parameters.spikeDetection

    units.unit = {};
    
    for u = 1:numel(Session.parameters.units)
        for f = fieldnames(Session.parameters.units(u))',
            f = f{1};
            if isnumeric(Session.parameters.units(u).(f))
                units.unit{u}.(f).Text = num2str(Session.parameters.units(u).(f));
                units.unit{u}.(f).Attributes.type = 'Numeric';                                
            else
                units.unit{u}.(f).Text = Session.parameters.units(u).(f);
                units.unit{u}.(f).Attributes.type = 'Char';                                                
            end
        end
    end    

    par.parameters.units = units;

    copyfile(fullfile(Session.spath,[Session.name,'.xml']),...
             fullfile(Session.spath,[Session.name,'.xml.bkp']));
    try,        
        struct2xml(par,fullfile(Session.spath,[Session.name,'.xml']));        
    catch err
        disp(err)
        %warning('MTA:MTASession:par2xml:AbortSave');
        %warning('MTA:MTASession:par2xml:RevertingToPreviousXML');
        copyfile(fullfile(Session.spath,[Session.name,'.xml.bkp']),...
                 fullfile(Session.spath,[Session.name,'.xml']));
    end
    
end
