function [distance,interact] = get_frootloop_distance(mouse_bodyparts,juv_bodyparts,varargin)
dist_thresh = 4/10*80; % 4 cm for miniscope, 10 cm per 80 pixels
mintime = 10; % frames
ignorepts = 5; % ignore small interaptions
juv_var = {'frootloop'};
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'dist_thr')
        dist_thresh = varargin{v+1};
    elseif  strcmpi(varargin{v},'mintime')
        mintime = varargin{v+1};
    elseif  strcmpi(varargin{v},'ignorepts')
        ignorepts = varargin{v+1};
    elseif  strcmpi(varargin{v},'discard_frames')
        discard_frames = varargin{v+1};
    elseif  strcmpi(varargin{v},'juv_var')
        juv_var = varargin{v+1};
    end
end
% get relative distance between test mouse and juvenile from bodyparts
distance = struct(); interact = struct();
mouse_var = {'implant','leftear','rightear'}; 
[mouse_position] =  get_bodyparts_position(mouse_bodyparts,mouse_var);
if ~isempty(discard_frames)
    mouse_position(discard_frames,:) = nan;
end
[frootloop_position] =  get_bodyparts_position(juv_bodyparts,juv_var);
distance.mouse_frootloop = arrayfun(@(x) pdist([mouse_position(x,:);frootloop_position(x,:)]),1:size(mouse_position,1));
interact.nose_frootloop = distance.mouse_frootloop;
%% nose(test mouse)- frootloop distance
% use implant position as test mouse nose
mouse_var = {'implant','leftear','rightear'}; 
[mouse_position] =  get_bodyparts_position(mouse_bodyparts,mouse_var);
if ~isempty(discard_frames)
    mouse_position(discard_frames,:) = nan;
end
[frootloop_position] =  get_bodyparts_position(juv_bodyparts,juv_var);
distance.nose_frootloop = arrayfun(@(x) pdist([mouse_position(x,:);frootloop_position(x,:)]),1:size(mouse_position,1));


%% butt (test mouse) - nose (juvenile) distance
mouse_var = {'tailbase'}; 
[mouse_position] =  get_bodyparts_position(mouse_bodyparts,mouse_var);
[frootloop_position] =  get_bodyparts_position(juv_bodyparts,juv_var);
distance.butt_frootloop = arrayfun(@(x) pdist([mouse_position(x,:);frootloop_position(x,:)]),1:size(mouse_position,1));

%% distance threshold: taken as the width of the juvenile's head
% leftear = get_bodyparts_position(mouse_bodyparts,{'leftear'});
% rightear = get_bodyparts_position(mouse_bodyparts,{'rightear'});
% dist_thresh = 2*nanmean(arrayfun(@(x) pdist([rightear(x,:);leftear(x,:)]),1:size(mouse_position,1)));
%% active investigating episodes
% 1. test mouse need to be close to juvenile
% 2. nose-nose distance should be shorter than butt-nose distance
init_trace = zeros(1,size(mouse_position,1));
interact.sniffs = init_trace;
interact.sniffs(distance.nose_frootloop<dist_thresh) = 1;
interact.proximal = init_trace; 
interact.proximal(distance.mouse_frootloop<dist_thresh) = 1;
%% approaching episodes
% 1. nose-nose and nose-butt distance should be decreasing
% 2. nose-nose distance should be shorter than butt-nose distance
temp = [diff(distance.nose_frootloop) 0];
% get episodes where the distance between the two mice continueously
% decreases
temp_approach = pt_continuousabove(-temp,0,0,mintime,Inf,ignorepts);
temp_approach_trace = init_trace;
for t = 1:size(temp_approach,1)
    temp_approach_trace(temp_approach(t,1):temp_approach(t,2)) = 1;
end
interact.approaching = init_trace;
interact.approaching(temp_approach_trace&distance.nose_frootloop<distance.butt_frootloop) = 1;

%% escaping episodes
% 1. nose-nose and nose-butt distance should be increasing
% 2. nose-nose distance should be longer than butt-nose distance
temp_escape = pt_continuousabove(temp,0,0,mintime,Inf,ignorepts);
temp_escape_trace = init_trace;
for t = 1:size(temp_escape,1)
    temp_escape_trace(temp_escape(t,1):temp_escape(t,2)) = 1;
end
interact.escaping = init_trace;
interact.escaping (temp_escape_trace&distance.nose_frootloop>distance.butt_frootloop&distance.nose_frootloop>dist_thresh) = 1;

%% filter out frames where the juvenile was not present
% [presence_frames,avg_lkl] = get_mouse_presence_frames(juv_bodyparts,0.7);
% presence_flag_trace = init_trace; presence_flag_trace(presence_frames) = 1;
% interact = structfun(@(x)x.*presence_flag_trace,interact,'UniformOutput',0);
interact = structfun(@(x)get_continuousabove_trace(x,0,ignorepts),interact,'UniformOutput',0);

% figure; hold on
% fds = fields(interact);
% arrayfun(@(x)plot(x+.7*interact.(fds{x}),'DisplayName',strrep(fds{x},'_',' ')),1:numel(fds))
% legend





end

%  dist_fun =@(x,y) sqrt((x(1)-y(1))^2+(x(2)-y(2))^2);
