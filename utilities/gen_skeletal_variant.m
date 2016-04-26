function [xyz] = gen_skeletal_variant(Trial,varargin)
[slist] = DefaultArgs(varargin,{'hand_labeled'},true);

slist = SessionList(slist);

local = false;
overwrite = true;

% Create Spline Spines
r1jid='';
for s = 1:numel(slist),
    Trial = MTATrial.validate(slist(s));
    file_preproc = fullfile(Trial.spath,[Trial.filebase '.pos.trb.t.mat']);
    if ~exist(file_preproc,'file')||overwrite,
        if ~local,
            system(['MatSubmitLRZ --config lrzc_hugemem.conf '...
                    r1jid ' -l ' Trial.name ...
                    ' transform_rigidBody']);
        else
            transform_rigidBody(slist(s));
        end
    end
end


