function [spont_trace] = get_spont_trace(input_trace,all_event_frames,peri_event_frames,min_spont_period)
% get peri-event frames to remove
exclud_frames = cell2mat(arrayfun(@(x)(x+[round(-peri_event_frames/2):1:round(peri_event_frames/2)]),all_event_frames,'un',0));
exclud_frames = unique(exclud_frames(:));
exclud_frames(exclud_frames>length(input_trace))=[]; exclud_frames(exclud_frames<=0)=[];
% remove spontaneous periods that are too short
spont_frames = ones(size(input_trace));
spont_frames (exclud_frames)=0;
long_spont_periods = pt_continuousabove(spont_frames,0,.5,min_spont_period,Inf,0);% smooth traces
% concatenate spontaneous periods
spont_frames = zeros(size(input_trace));
for i = 1:size(long_spont_periods)
    spont_frames(long_spont_periods(i,1):long_spont_periods(i,2)) = 1;
end
spont_trace = input_trace; spont_trace(spont_frames<1)= nan;
end

