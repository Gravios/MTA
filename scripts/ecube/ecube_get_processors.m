function xml = ecube_get_processors(xml)

sig = true([numel(xml.SETTINGS.SIGNALCHAIN),1]);
if numel(xml.SETTINGS.SIGNALCHAIN)==1, 
    xml.SETTINGS.SIGNALCHAIN = {xml.SETTINGS.SIGNALCHAIN};
else
    for s = 1:numel(xml.SETTINGS.SIGNALCHAIN),
        if numel(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR)==1,
            if ~isfield(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR,'CHANNEL'),
                sig(s) = false;
                continue;
            end
            xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR = {xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR};
        else
            proc = true([numel(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR),1]);
            for p = 1:numel(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR),
                if ~isfield(xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR{p},'CHANNEL'),
                    proc(p) = false;
                end
            end
            xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR = xml.SETTINGS.SIGNALCHAIN{s}.PROCESSOR(proc);
        end
    end
    xml.SETTINGS.SIGNALCHAIN = xml.SETTINGS.SIGNALCHAIN(sig);
end
