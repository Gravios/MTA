function generate_session_report(Session)


global MTA_PROJECT_PATH

hfig = figure(20161109)

document.format = 'A4'

document.name = [Session.filebase,'_rpt.pdf'];
document.path = fullfile(MTA_PROJECT_PATH,Session.name)
                 

                 
                 
switch document_size
  case 'A4'
    document.height = 297;
    document.width = 210;
end

document.orientation = 'portrait';

% page1 ---------------------------------------------------------------------



% Session Meta Data

% XYZ occupancy summary 


% page2 ---------------------------------------------------------------------