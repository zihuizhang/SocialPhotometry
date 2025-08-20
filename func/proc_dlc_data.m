function [mouse_position,chamber_frames,chamber_masks,cup_centroids,...
    mouse_bodyparts,juv_bodyparts,juv_presence_frames] = proc_dlc_data(this_file,this_session_type,frame_size)
%% get mouse position, chamber masks, and cup positions from DLC result file
% ZZ 2022
[dlc_vars,chamber_names] = set_dlc_vars(this_session_type);
juv_presence_frames = [];
%% read dlc output csv
thresh = 0.9;
if ischar(this_file)
data = csvread(this_file,3,1);
[~,this_name] = fileparts(this_file);

elseif isnumeric(this_file)
    data = this_file;
    this_name = '';
end
num_vars = numel(dlc_vars); num_frames = size(data,1);
% check data size. each var should have x,y,likelihood
if size(data,2)~= num_vars*3
    disp('dlc vars and data mismatch!!')
end
% interpolate points with low likelihood
for v = 1:num_vars
    this_lkhd = data(:,3*v);
    this_x = data(:,3*(v-1)+1);
    this_y = data(:,3*(v-1)+2);
    this_bad_frames = find(this_lkhd<thresh);
    if numel(this_bad_frames)>length(this_lkhd)*0.5
        continue
    end
    this_x(this_bad_frames) = nan;
    this_y(this_bad_frames) = nan;
    data(:,3*(v-1)+1) = interp1_nans(this_x);
    data(:,3*(v-1)+2) = interp1_nans(this_y);
end
%% get mouse position
% average ear, tail and implant positions
if contains(this_session_type,'HomecageTwoMice')||contains(this_session_type,'SIMBA_Aggression')||contains(this_session_type,'HomecageIntruder')||contains(this_session_type,'JuvenileInteractionBox61')
    mouse_vars = {'Nose_1','Ear_left_1','Ear_right_1','Tail_base_1','Center_1'};
else
    mouse_vars = {'leftear','rightear','tailbase','implant'};
end
[mouse_position,mouse_bodyparts,mouse_heatmap] =  get_mouse_position(data, dlc_vars,mouse_vars,frame_size);

%% get juvenile position (for 'TwoMice' data)
if contains(this_session_type,'openfieldTwoMice')
    juv_vars =  {'juv_snout','juv_leftear','juv_rightear','juv_tailbase'};
    [juv_position,juv_bodyparts,juv_heatmap] =  get_mouse_position(data, dlc_vars,juv_vars,frame_size);
    juv_presence_frames = get_mouse_presence_frames(juv_bodyparts,0.6);
elseif contains(this_session_type,'HomecageTwoMice')||contains(this_session_type,'SIMBA_Aggression')||contains(this_session_type,'JuvenileInteractionBox61')
    juv_vars =  {'Nose_2','Ear_left_2','Ear_right_2','Tail_base_2','Center_2'};
    [juv_position,juv_bodyparts,juv_heatmap] =  get_mouse_position(data, dlc_vars,juv_vars,frame_size);
%     juv_presence_frames = get_mouse_presence_frames(juv_bodyparts,0.6);
elseif contains(this_session_type,'homecageFrootloop')||contains(this_session_type,'openfieldFrootloop')
    juv_vars =  {'frootloop'};
    [juv_position,juv_bodyparts,juv_heatmap] =  get_mouse_position(data, dlc_vars,juv_vars,frame_size);
    juv_presence_frames = get_mouse_presence_frames(juv_bodyparts,0.6);
else
    juv_position = []; juv_bodyparts = []; juv_heatmap = []; juv_presence_frames = [];
end

%% get chamber edges
chamber_edge_names = {'lowerleft','upperleft','upperright','lowerright',};
chamber_masks = struct();

% get the average position from all frames
for c = 1:numel(chamber_names)
    this_chamber = chamber_names{c};
    this_vars = cellfun(@(x)[x,this_chamber],chamber_edge_names,'un',0);
    [~,~,var_idx]=intersect(this_vars,dlc_vars,'stable');
    this_x = nanmean(data(:,3.*(var_idx-1)+2),1);
    this_y = nanmean(data(:,3.*(var_idx-1)+1),1);
    this_mask =  poly2mask(this_x,this_y,frame_size(1),frame_size(2));
    chamber_masks.(['chamber' this_chamber]) = this_mask;
end

% add the center chamber mask for threechamber
if strcmp(this_session_type,'threechamber')
    this_chamber = 'C';
    [~,~,var_idx]=intersect({'upperleftB','lowerleftA','lowerrightA','upperrightB'},dlc_vars,'stable');
    this_x = nanmean(data(:,3.*(var_idx-1)+2),1);
    this_y = nanmean(data(:,3.*(var_idx-1)+1),1);
    this_mask =  poly2mask(this_x,this_y,frame_size(1),frame_size(2));
    chamber_masks.(['chamber' this_chamber]) = this_mask;
end

%% get cup centroids (if any)
cup_centroid_names= {'cupA','cupB',};
cup_centroids = struct();
for c = 1:numel(cup_centroid_names)
    this_cup = cup_centroid_names{c};
    [~,~,var_idx]=intersect(this_cup,dlc_vars,'stable');
    if~isempty(var_idx)
        this_x = nanmean(data(:,3.*(var_idx-1)+2),1);
        this_y = nanmean(data(:,3.*(var_idx-1)+1),1);
        cup_centroids.(this_cup) = [this_x this_y];
    end
end
%% get frames spent in each chamber
chamber_frames = struct();
fds = fields(chamber_masks);
for f = 1:numel(fds)
    fd = fds{f};

    % interset indices
    [this_mask_x,this_mask_y] = find(chamber_masks.(fd)>0);
    this_event_trace = zeros(1,num_frames);
%     this_posisions = arrayfun(@(x)sub2ind(frame_size,mouse_position(x,2),mouse_position(x,1)),1:size(mouse_position,1));
%     this_event_trace = arrayfun(@(x)any(this_mask(:) == x),this_posisions(:));

    this_event_trace = arrayfun(@(x,y)(~isempty(intersect(this_mask_x,x)))&(~isempty(intersect(this_mask_y,y))),mouse_position(:,2),mouse_position(:,1),'UniformOutput',1);
    chamber_frames.(fd) = this_event_trace;
end
%% plot
figure('name','dlc result','units','normalized','outerposition',[0 0 1 1])
fds = fields(chamber_masks);
num_plot_cols = numel(fds)+1;
for f = 1:numel(fds)
    fd = fds{f};
    this_mask = chamber_masks.(fd);
    subplot(1,num_plot_cols,f); hold on

    if contains(this_session_type,'HomecageTwoMice')|| contains(this_session_type,'SIMBA_Aggression')||...
            contains(this_session_type,'openfieldTwoMice')||contains(this_session_type,'JuvenileInteractionBox61')  % session types where the chamber edges were not labelled
        plot(mouse_position(:,1),mouse_position(:,2),'color',[.7 .7 .7])
    else
        imagesc(this_mask); title(fd); 
        plot(mouse_position(chamber_frames.(fd),1),mouse_position(chamber_frames.(fd) ,2),'color',[.5 .5 .5])

    end

    if (~isempty(juv_position))&&(~numel(juv_presence_frames)>2*size(juv_position,1)/3) % otherwise not a juvenile session
            plot(juv_position(juv_presence_frames,1),juv_position(juv_presence_frames ,2),'color','y','LineStyle',':')
    end

    xlim([0 frame_size(1)]);ylim([0 frame_size(1)])
    set(gca,'YDir','reverse')
    axis square
end
subplot(1,num_plot_cols,f+1); hold on
imagesc(mat2gray(mouse_heatmap)); colormap(gray)
% plot(mouse_position(:,1),mouse_position(:,2),'w')
try
    scatter(cup_centroids.cupA(1),cup_centroids.cupA(2),'markerfacecolor','r','markeredgecolor','w');
    scatter(cup_centroids.cupB(1),cup_centroids.cupB(2),'markerfacecolor','b','markeredgecolor','w');
end
xlim([0 frame_size(1)]);ylim([0 frame_size(1)])
axis square
set(gca,'YDir','reverse')
title('mouse and cup(s)')
suptitle(strrep(this_name,'_',' '))


%% 

end

