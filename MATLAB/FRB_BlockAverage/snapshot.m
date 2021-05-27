% function description : get candidate snapshot from the raw data

% raw_data : The raw data that we cut snapshot from
% target_length: The time length of target signal segment (unit: s)
% raw_length: length of raw signal
% extra_factor: A factor to save a snapshot longer than target signal length (to add some buffer so as not to miss any target signal)
% trigger_time: The time of triggering (unit: s)
% output: candidate snapshot data (time domain)
% fs: sampling frequency of raw data
% snapshot_flag: flag of snapshot


function output = snapshot(raw_data,target_length,raw_length,extra_factor,trigger_time,fs,snapshot_flag)

% get start time of snapshot
snapshot_start = trigger_time - target_length - target_length*extra_factor;
start_index = ceil(snapshot_start/(1/fs));               % transfer start time to vector index
if snapshot_start <= 0                                   
    snapshot_start = 0;
    start_index = 1;
end

% get end time of snapshot
snapshot_end = trigger_time + target_length*extra_factor;
end_index = ceil(snapshot_end/(1/fs));                    % transfer end time to vector index
if snapshot_end > raw_length                               % if snapshot_end exceed the length of the raw signal
    snapshot_end = raw_length;
    end_index = snapshot_end/(1/fs);
end

% get snapshot
output = raw_data(start_index:end_index);
snapshot_t = 0:1/fs:(end_index-start_index)*1/fs;
plotPlus(output,snapshot_t,fs,'Candidate snapshot');

end



