

% INIT Feature
fet = MTADfet(Trial.spath,                     ... Path              (String)
              [],                              ... File Name         (String)
              [],                              ... Data              (Numeric)
              Trial.xyz.sampleRate,            ... Sample Rate       (Numeric)
              Trial.sync.copy,                 ... Sync Periods      (MTADepoch)
              Trial.sync.data(1),              ... Sync Origin       (Numeric)
              [],                              ... Model             (MTAModel)
              'TimeSeries',                    ... DataType          (String)
              [],                              ... File Extension    (String)
              'body_referenced_markers',       ... Name              (String)
              'bref',                          ... Label             (String)
              'b'                              ... Key               (String)
);


% PREPROC xyz
xyz = preproc_xyz(Trial,procOpts0;
