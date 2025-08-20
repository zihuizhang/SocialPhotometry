function [mouse_heatmap] = get_mouse_heatmap(mouse_position,frame_size)

mouse_heatmap = zeros(frame_size);
for p = 1:size(mouse_position,1)
    if sum(isnan(mouse_position(p,:)))<1&&(mouse_position(p,2)<=frame_size(1))&&(mouse_position(p,1)<=frame_size(2)&&min(mouse_position(p,:))>0)
        mouse_heatmap(mouse_position(p,2),mouse_position(p,1)) =1+mouse_heatmap(mouse_position(p,2),mouse_position(p,1));
    else
        continue
    end
end
mouse_heatmap = imgaussfilt(mouse_heatmap,8); % blur heat map

end