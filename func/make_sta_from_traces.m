function [ mean_dff_sta_traces,dff_sta,mean_sta_traces,raw_sta_traces,stim_frames,df_sta,mean_df_sta] = make_sta_from_traces( traces,frames_with_stim,pre_frames,post_frames,baseline_frames )
%
if(size(traces,2)<size(traces,1))
    traces = traces';
end

mean_dff_sta_traces = [];
dff_sta = [];
mean_sta_traces = [];
df_sta = [];
mean_df_sta = [];
raw_sta_traces = [];
stim_frames = [];


if (isempty(frames_with_stim))
    warning('no stim found, make_sta_from_traces return')
    return 
end

trace_length = size(traces,2);
stim_frames = frames_with_stim((frames_with_stim<trace_length-post_frames));
raw_sta_traces = make_sta_traces(traces, stim_frames, pre_frames, post_frames);
dff_sta = nan(size(raw_sta_traces));
df_sta = dff_sta;

% get baselined sta traces
sta_baseline = nanmean (raw_sta_traces(:,baseline_frames),2);
for m = 1:size(sta_baseline,1) 
        dff_sta(m,:) = (raw_sta_traces(m,:)-sta_baseline(m))./sta_baseline(m);
        df_sta(m,:) = (raw_sta_traces(m,:)-sta_baseline(m)); 
end

mean_dff_sta_traces = squeeze(nanmean(dff_sta, 1));  % targets trace mean, averaged over stims
mean_sta_traces = squeeze(nanmean(raw_sta_traces, 1));
mean_df_sta = squeeze(nanmean(df_sta, 1));
dff_sta = squeeze(dff_sta);
df_sta = squeeze(df_sta);
% median_dff_sta_traces = squeeze(nanmedian(dff_sta, 2));


end

