function processorList = ecube_select_ephys_sources(xml)

processorList = {};
sig = true([numel(xml.SETTINGS.SIGNALCHAIN),1]);
if numel(xml.SETTINGS.SIGNALCHAIN)==1, 
    xml.SETTINGS.SIGNALCHAIN = {xml.SETTINGS.SIGNALCHAIN};
else
    for s = 1:numel(xml.SETTINGS.SIGNALCHAIN),
        if numel(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR)==1,
            if isfield(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR,'CHANNEL'),
                processorList{end+1} = xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR.Attibutes.NodeId;
            end
        else
            for p = 1:numel(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR),
            if isfield(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR{p},'CHANNEL'),
                processorList{end+1} = xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR{p}.Attibutes.NodeId;
            end
        end
    end
end
