function [bodyparts_position] =  get_bodyparts_position(mouse_bodyparts,bodyparts)
% get average position of specific body parts from mouse_bodayparts
% structure saved out by 'get_mouse_position'
data = cell2mat(struct2cell(mouse_bodyparts)');
xidx =cell2mat(cellfun(@(x)find(contains(fields(mouse_bodyparts),[x,'_x'])),bodyparts,'UniformOutput',0));
this_mouse_x = mean(data(:,xidx),2);
yidx =cell2mat(cellfun(@(x)find(contains(fields(mouse_bodyparts),[x,'_y'])),bodyparts,'UniformOutput',0));
this_mouse_y =  mean(data(:,yidx),2);
bodyparts_position = round([this_mouse_x this_mouse_y]);
bodyparts_position(bodyparts_position==0) = 1; bodyparts_position(bodyparts_position>512) = 512;
end