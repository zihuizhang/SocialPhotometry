function [opt] = set_opt0(session_data,varargin)
sta_pre_sec = 20;
sta_post_sec = 60;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'sta_pre_sec')
        sta_pre_sec = varargin{v+1};
    end
    if strcmpi(varargin{v},'sta_post_sec')
        sta_post_sec = varargin{v+1};
    end
end
% default parameters for single session
 opt = struct();
    opt.red_channel_name = 'Data_stream_560';
    opt.green_channel_name = 'Data_stream_465';
    opt.uv_channel_name = 'Data_stream_405';
    opt.ds_factor = 3; % added 20220414. take moving average of TDT recording to speed up
    opt.frame_rate = round(session_data.sampling_frequency/opt.ds_factor);
    opt.discard_first_frames = round(opt.frame_rate*5);
    opt.sta_pre_sec = sta_pre_sec;
    opt.sta_post_sec = sta_post_sec;
    opt.baseline_win_start_sec = 8; % sec before stim
    opt.baseline_duration_sec = 4;
    opt.pre_frames = round(opt.sta_pre_sec*opt.frame_rate);
    opt.post_frames = round(opt.sta_post_sec*opt.frame_rate);
    opt.baseline_frames = opt.pre_frames - round(opt.baseline_win_start_sec*opt.frame_rate) +[1:round(opt.frame_rate*opt.baseline_duration_sec)];
    opt.num_frames = round(length(downsample(session_data.(opt.red_channel_name),opt.ds_factor)));
    opt.session_duration = opt.num_frames/opt.frame_rate;
    opt.filt_param = 0.01; % for lowess smoothing (before taking correlation)
    opt.sta_corr_range_sec = [-2 2];
    opt.sta_corr_range = round(opt.sta_corr_range_sec.*opt.frame_rate+opt.pre_frames);
    % behav parameters
    opt.behav_frame_rate = round(opt.frame_rate); % Hz
    opt.behav_pre_sec= 10;
    opt.behav_post_sec = 10;
    opt.behav_pre_frames = opt.behav_frame_rate*opt.behav_pre_sec;
    opt.behav_post_frames = opt.behav_frame_rate*opt.behav_post_sec;
    opt.behav_xticks = [-opt.behav_pre_frames:1:0:opt.behav_post_frames]./opt.behav_frame_rate;
    opt.frustration_start = 300; % sec after first reward
    opt.behav_num_frames = ceil(opt.num_frames/opt.frame_rate*opt.behav_frame_rate);
    opt.event_binsize_sec = 30;
    opt.plot_sta_fds = {'uv','green',};
    opt.plot_sta_names = {'iso','GRAB'};
    opt.location = 'NAc';
    opt.session_type = '';

end

