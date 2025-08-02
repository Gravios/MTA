function generate_session_report_page1(Session)

hfig = figure(20161109)

document_size = 'A4'

document = struct('name', ,...
                  'path', ,...
switch document_size
  case 'A4'
    document.height = 297