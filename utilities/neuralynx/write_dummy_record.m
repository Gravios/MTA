function write_dummy_record(fid,timestamp,chanId,sampleRate,nSamples)
    fwrite(fid, uint64(timestamp), 'uint64');     % block timestamp
    fwrite(fid, uint32(chanId),    'uint32');     % block channel
    fwrite(fid, uint32(sampleRate),'uint32');     % block SAMPLE_RATE
    fwrite(fid, uint32(nSamples),  'uint32');     % block samples
    fwrite(fid, uint16(zeros([1,512])),'uint16'); % block dummy data
end
