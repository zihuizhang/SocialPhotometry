function [this_attack_start_times,this_attack_stop_times] = filt_event_timestamps(this_attack_start_times,this_attack_stop_times,min_dur, max_isi)
    % discard events that are too short and merge events close together
    % merge close events
    isi = this_attack_start_times(2:end) - this_attack_stop_times(1:end-1);
    close_event_idx = find(isi<max_isi);
    if ~isempty(close_event_idx)
        this_attack_start_times(close_event_idx+1) = [];
        this_attack_stop_times(close_event_idx) = [];
    end
    
    % discard short events
    dur = this_attack_stop_times-this_attack_start_times;
    this_attack_stop_times(dur<min_dur) = [];
    this_attack_start_times(dur<min_dur) = [];
    

end

