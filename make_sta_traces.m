function [sta_traces,df_sta_traces,dff_sta_traces] = make_sta_traces(input_traces, stim_frames, pre_frames, post_frames,varargin)
baseline_frames = 1:pre_frames-1; % use the entire pre stim window for baselining unless specified
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'baseline_frames')
        baseline_frames = varargin{v+1};
    end
end

stim_frames((stim_frames-pre_frames)<1) = [];
stim_frames((stim_frames+post_frames)>length(input_traces)) = [];
sta_traces = zeros(length(stim_frames), pre_frames+post_frames+1);
df_sta_traces = sta_traces;
dff_sta_traces = sta_traces;

for stim = 1:length(stim_frames)
    try
        sta_traces( stim, :) = input_traces( stim_frames(stim)-pre_frames:stim_frames(stim)+post_frames);
        df_sta_traces( stim, :) = sta_traces( stim, :)- mean( sta_traces( stim,baseline_frames));
        dff_sta_traces( stim, :) = df_sta_traces( stim, :)./mean( sta_traces( stim, baseline_frames));
    catch
        sta_traces( stim, :) = NaN;
    end
end

