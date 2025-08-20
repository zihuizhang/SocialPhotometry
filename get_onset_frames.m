function [attack_onset_frames,attack_offset_frames] = get_onset_frames(attack_frames)
if isempty(attack_frames)
    attack_onset_frames = []; attack_offset_frames = [];
else
    onset_frame_idx = find(diff(attack_frames)>1)+1;
    onset_frame_idx(onset_frame_idx>length(attack_frames)) = [];
    attack_onset_frames = [attack_frames(1);attack_frames(onset_frame_idx)-1];
    
    offset_frame_idx = diff(attack_frames)>1;
    attack_offset_frames = [attack_frames(offset_frame_idx);attack_frames(end)];
end
end

