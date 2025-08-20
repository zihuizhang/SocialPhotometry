function [presence_frames,avg_lkl] = get_mouse_presence_frames(juv_bodyparts,thresh)
% infer when the mouse is present in the video from the avg. bodypart
% likelihood. Mouse is probably not present if likelihood<thresh
fds = fields(juv_bodyparts);
fds = fds(contains(fds,'_lkl'));
lkl_vals = cellfun(@(x)juv_bodyparts.(x)(:),fds,'UniformOutput',false);
lkl_vals = cell2mat(lkl_vals');
avg_lkl = mean(lkl_vals,2);
episodes = pt_continuousabove(avg_lkl,0,thresh,60,Inf,20);
presence_frames = [];
for t = 1:size(episodes,1)
    presence_frames = [presence_frames episodes(t,1):episodes(t,2)];
end
%%
%  figure; plot(avg_lkl)
end