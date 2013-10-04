function linkSession(session_name,target_path,paths)
   if exist([target_path.xyz session_name],'file')&&exist([target_path.nlx session_name],'file'),

       mkdir([ paths.xyz session_name ]);
       mkdir([ paths.nlx session_name ]);
       mkdir([ paths.analysis session_name ]);

       subsessions = dir([target_path.xyz session_name]);
       subsessions = subsessions(3:end);
       for i = 1:length(subsessions),
           mkdir([ paths.xyz session_name '/' subsessions(i).name]);
           session_parts = dir([target_path.xyz session_name '/' subsessions(i).name]);
           session_parts = session_parts(3:end);
           for j = 1:length(session_parts)
               system(['ln -s ' target_path.xyz session_name '/' subsessions(i).name '/' session_parts(j).name ' ' paths.xyz session_name '/' subsessions(i).name '/' session_parts(j).name]);
           end
       end

       nlx_file = dir([target_path.nlx session_name]);
       nlx_file = nlx_file(3:end);
       for j = 1:length(nlx_file)
           system(['ln -s ' target_path.nlx session_name '/' nlx_file(j).name ' ' paths.nlx session_name '/' nlx_file(j).name]);
       end

       analysis_file = dir([target_path.analysis session_name]);
       analysis_file = analysis_file(3:end);
       for j = 1:length(analysis_file)
           system(['ln -s ' target_path.analysis session_name '/' analysis_file(j).name ' ' paths.analysis session_name '/' analysis_file(j).name]);
       end

  end
end
