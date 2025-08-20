function [mouse_position,mouse_bodyparts,mouse_heatmap] =  get_mouse_position(data, dlc_vars,mouse_vars,frame_size, varargin)
mouse_vars_savenames = mouse_vars;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'mouse_vars_names')
        mouse_vars_savenames = varargin{v+1};
    end
end
[~,var_idx]=intersect(dlc_vars,mouse_vars);
this_mouse_x = mean(data(:,3.*(var_idx-1)+2),2);
this_mouse_y = mean(data(:,3.*(var_idx-1)+1),2);
mouse_position = round([this_mouse_x this_mouse_y]);
mouse_position(mouse_position==0) = 1; mouse_position(mouse_position>512) = 512;
[mouse_heatmap] = get_mouse_heatmap(mouse_position,frame_size);
mouse_bodyparts = struct();
% save all body part positions
for v = 1:numel(mouse_vars)
    vv = mouse_vars{v}; [~,var_idx]=intersect(dlc_vars,vv);
    vv_save = mouse_vars_savenames{v};
    mouse_bodyparts.([vv_save '_x']) = data(:,3.*(var_idx-1)+2);
    mouse_bodyparts.([vv_save '_y']) = data(:,3.*(var_idx-1)+1);
    mouse_bodyparts.([vv_save '_lkl']) = data(:,3.*(var_idx-1)+3);

end

end