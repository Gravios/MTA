
mode = 'build_TimeStampMap';
fbase = 'jg05-20120312-01-cof-nlx-0';

ncs_ext = '.ncs';
ts_ext = '.timestamps';
dat_ext = '.dat';

   
switch mode

  case 'extract_Data&TimeStamps'

for f = 1:96,
    if ~(f > 9),
        fnum = ['0',num2str(f)];
    else
        fnum = num2str(f);
    end

    fid = fopen([fbase,fnum,ncs_ext],'r');
    oid = fopen([fbase,fnum,ts_ext],'a+');
    did = fopen([fbase,fnum,dat_ext],'a+');

    fseek(fid,0,'eof');
    eofpos = ftell(fid);

    fseek(fid,0,'bof');
    fseek(oid,0,'bof');
    fseek(did,0,'bof');

    hdr = []; freq=[]; rsize=[];

    % NCS file header
    hdr = fread(fid,16*1024,'char=>char');
    status = fseek(fid,16*1024,'bof');

    i=1;
    ts = fread(fid,1,'uint64=>uint64');
    fwrite(oid,ts,'integer*8');
    status = fseek(fid,4,'cof'); 
    freq = single(fread(fid,1,'uint32=>uint32'));
    rsize = single(fread(fid,1,'uint32=>uint32'));
    dat = fread(fid,rsize,'int16=>int16');
    fwrite(did,dat,'integer*2');
    %fseek(fid,1024,'cof'); % skip record samples
    tic
    while rsize==512&&freq==32556&&ftell(fid)~=eofpos
        i=i+1;
        ts = fread(fid,1,'uint64=>uint64');
        fwrite(oid,ts,'integer*8');
        fseek(fid,4,'cof');
        freq = single(fread(fid,1,'uint32=>uint32')); 
        rsize = single(fread(fid,1,'uint32=>uint32'));
        dat = fread(fid,rsize,'int16=>int16');
        fwrite(did,dat,'integer*2');
        %fseek(fid,1024,'cof'); % skip record samples
    end
    toc
    fclose(fid)
    fclose(oid)
    fclose(did)
end


case 'build_TimeStampMap';

for f = 1:96,
    if ~(f > 9),
        fnum = ['0',num2str(f)];
    else
        fnum = num2str(f);
    end
    oid(f) = fopen([fbase,fnum,ts_ext],'a+');
    fseek(oid(f),0,'eof');
    oid_eof = ftell(oid(f));
    fseek(oid(f),0,'bof');
end

tss = [];
i = 1;
for f = 1:96,
tss(i,f) = fread(oid(f),1,'uint64=>uint64');
end
unq_ts = unique(tss(i,:));
if numel(unq_ts)>1,
    error('timestamps misalignment at start')
    %min_ts = min(unq_ts)
end

%while ftell(oid)~=oid_eof,

i=i+1;
%end



function check_future_ts(fid,)

end


for f = 1:96,
    fclose(oid(f));
end

did = fopen('jg05-20120312-01-cof-nlx-025.dat','r');
fseek(did,0,'eof');
did_eof = ftell(oid);
fseek(did,0,'bof');
data=[];
i = 1;
while ftell(oid)~=did_eof,
data((1+(i-1)*512):((i-1)*512+512)) = fread(did,512,'int16=>int16');
i=i+1;
end
fclose(did)