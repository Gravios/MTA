% req20210325
%    Tags: cabels.gl visualization
%    Status: Active
%    Type: Diagnostic
%    Author: Justin Graboski
%    Final_Forms: NA
%    Project: General
%    Description: Generate json file with xyz data for all markers for input file to cables.gl


Trial = MTATrial.validate('Ed03-20140625.cof.all');

xyz = Trial.load('xyz');

fid = fopen('/storage/gravio/test6.json','w');

fprintf(fid,'[\n');
for t = 1:5000,
fprintf(fid,['\t{"X":[%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f],', ...
             ' "Y":[%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f],', ...
             ' "Z":[%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f]},\n'], ...
             reshape(xyz(t,:,:),[],1));
end

fprintf(fid,['\t{"X":[%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f],', ...
             ' "Y":[%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f],', ...
             ' "Z":[%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,%4.2f]}\n'], ...
             reshape(sq(xyz(t+1,:,:)),[],1));
fprintf(fid,']');


fclose(fid);


fid = fopen('/storage/gravio/test7.json','w');
fprintf(fid,'[\n');
for t = 1:5000,
fprintf(fid,['\t{"MarkerXYZ":[%4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, ', ...
             '%4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f ,', ...
             '%4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f]},\n'], ...
             reshape(sq(xyz(t,:,:))',[],1));
end
fprintf(fid,['\t{"MarkerXYZ":[%4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, ', ...
             '%4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f ,', ...
             '%4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f, %4.2f,%4.2f,%4.2f]}\n'], ...
             reshape(sq(xyz(t+1,:,:))',[],1));
fprintf(fid,']');
fclose(fid);