% simple pipeline for photometry social interactions
% requires converted TDT mat files and DLC output spreadsheets
% - pre-process .mat files converted from tdt files
% - get social interaction timestamps from DLC files (distance<5cm)
% - get interaction-aligned photometry traces
% ZZ 2024
% >>>> PLEASE KEEP THE NAMING FORMATS CONSISTENT <<<<<
%% set directories - CHANGE THE DIRECTORIES TO YOUR OWN!!
% matlab paths 
matlab_path = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry'; % CHANGE THIS
addpath(genpath(matlab_path)); cd(matlab_path)

% analysis plan
excel_file = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\2024_NAcInputs_Juvenile.xlsx';  % CHANGE THIS

% data paths
matfile_dir = 'E:\ZoeVideos\TDTMATFiles\juvenile'; % put TDT converted (.mat) files here  % CHANGE THIS
dlc_dir = 'E:\ZoeVideos\DLCprocdata\processed'; % put DLC output (.csv) files here  % CHANGE THIS
tdt_mat_path = 'E:\ZoeVideos\TDTProcData';  % processed data will be saved here % CHANGE THIS
mat_fileList = dir(fullfile(matfile_dir,'/*.mat'));
dlcfileList = dir(fullfile(dlc_dir));


%% select session type
DLC_model = 'JuvenileInteractionBox61'; stimulus_name = 'Social';
sheet_name = 'mPFC';
%% params
video_frame_rate = 20; % Hz; resolution for behav traces
min_social_dur = 2; % sec used 2 sec for GRAB5HT data
max_social_isi = 1; % sec
colors = set_colors();
Tprocess = readtable(excel_file,'Sheet',sheet_name);
IF_LOAD_PROC = 0;
pixelspermm = 1.7; % for box6.1
dist_thresh = pixelspermm*50; % in pixels, 5cm
event_types = cellfun(@(x)[stimulus_name x],{'Onset','Offset','Entry'},'un',0);
%% global params for photometry
global IF_GET_CORRELATION IF_LOAD_PROCTRACE IF_PLOT_RAW IF_ADD_FIELDS
IF_ADD_FIELDS = false; % load and add field to existing data (temporal use)
IF_GET_CORRELATION = false; % correlogram of red and green
IF_LOAD_PROCTRACE = true;
IF_PLOT_RAW = false;
%% loop through animals to process individual sessions
all_mat_files = {mat_fileList.name}';
num_sessions = size(Tprocess,1);
for a = 1:num_sessions
    this_animal = Tprocess.Mouse{a};
    this_animal_index_str = strrep(strrep(this_animal,'ZZ',''),'PJA','');
    this_animal_index_pad =  num2str(str2double(this_animal_index_str),'%03.f');
    this_date = Tprocess.Date(a,:);
    this_program = Tprocess.Program{a};
    this_TDTend = Tprocess.TDTend(a,:);
    save_struct = struct();
    this_start_time = Tprocess.SessionStartFrame(a,:)/video_frame_rate; % when the juvenile was introduced
    try
        this_TDTduration = Tprocess.SessionDuration(a,:); % sec
    catch
        this_TDTduration = [];
    end
    
    %% find TDT file
    this_mat_date = datestr(this_date,'yymmdd');
    if ~isnan(this_TDTend)
        this_mat_date = [this_mat_date '-' num2str(this_TDTend,'%06.f')];
    end
    if ~strcmp(this_animal_index_pad,'NaN')
        this_mat_animal = strrep(this_animal,this_animal_index_str,this_animal_index_pad);
    else
        this_mat_animal = this_animal;
    end
    this_mat_file = all_mat_files{cellfun(@(x)(contains(x,strrep(this_mat_animal,'ZZ','ZZH24-'))...
        ||contains(x,strrep(this_mat_animal,'ZZ','ZZH23-'))||contains(x,strrep(this_mat_animal,'ZZ','ZZH22-'))||contains(x,strrep(this_mat_animal,'PJA','PJA24-'))...
        )&&contains(x,this_mat_date),all_mat_files)};
    
    this_save_date = datestr(this_date,'yyyymmdd');
    
    %% load dlc to get aggression episodes
    % read dlc tables
    input_excel_name = dlcfileList(cellfun(@(x)contains(x,this_animal)&&(contains(x,datestr(Tprocess.Date(a,:),'yyyymmdd'))||contains(x,datestr(Tprocess.Date(a),'mmddyyyy')))...
        &&contains(x,this_program),{dlcfileList(:).name}));
    input_excel_path = [dlc_dir filesep input_excel_name.name];
    if isempty(input_excel_name)
        disp([datestr(Tprocess.Date(a,:),'yyyymmdd') ' ' this_animal ' ' this_program ' DLC csv not found!' ])
    end
    
    Tdata =readtable(input_excel_path);
    %% get locomotion and distance between mice
    [dlc_vars] = set_dlc_vars(DLC_model);
    [mouse_position,~,~,~,...
        mouse_bodyparts,juv_bodyparts] = proc_dlc_data(input_excel_path,DLC_model,[512,512]);
    close
    if strcmp(DLC_model, 'JuvenileInteractionBox61')
        [distance,interact] = get_intruder_distance(mouse_bodyparts,juv_bodyparts,'dist_thr',dist_thresh,'discard_frames',1:this_start_time*video_frame_rate);
    elseif  strcmp(DLC_model,'FrootloopBox61')
        [distance,interact] = get_frootloop_distance(mouse_bodyparts,juv_bodyparts,'juv_var','frootloop',...
            'dist_thr',dist_thresh,'discard_frames',1:this_start_time*video_frame_rate);
    elseif strcmp(DLC_model,'ToymouseBox61')
        [distance,interact] = get_frootloop_distance(mouse_bodyparts,juv_bodyparts,'juv_var','toymouse',...
            'dist_thr',dist_thresh,'discard_frames',1:this_start_time*video_frame_rate);
    end
    mouse_speed = arrayfun(@(x)pdist([mouse_position(x,:);mouse_position(x+1,:)],'euclidean'),1:size(mouse_position,1)-1)*video_frame_rate;
    mouse_distance = distance.mouse_juv;
    mouse_speed = mouse_speed./pixelspermm./10;% cm/s
    mouse_distance = mouse_distance./pixelspermm./10;% cm
    
    [event_onset_frames,event_offset_frames] = structfun(@(x)get_onset_frames(find(x>0)'),interact,'UniformOutput',false);
    
    % discard events before juvenille put in
    mouse_distance(1:this_start_time*video_frame_rate) = nan;
    
    % discard events that are too short and merge events close together
    [event_onset_frames.proximal,event_offset_frames.proximal] = filt_event_timestamps(event_onset_frames.proximal,event_offset_frames.proximal,...
        min_social_dur*video_frame_rate, max_social_isi*video_frame_rate);
    
    % save event timestamps to struct
    save_struct.([stimulus_name 'Onset']) = event_onset_frames.proximal./video_frame_rate;
    save_struct.([stimulus_name 'Offset']) = event_offset_frames.proximal./video_frame_rate;
    save_struct.([stimulus_name 'Entry']) = this_start_time;
    
    behav_timestamps = save_struct;
    save_struct.speed = mouse_speed;% cm/s
    save_struct.distance = mouse_distance;% cm
    
    %% process photometry data
    this_load_path = [tdt_mat_path filesep this_animal filesep this_save_date filesep this_program];
    if IF_LOAD_PROC&&exist([this_load_path filesep 'sta_struct.mat'],'file')
        sta_struct = load([this_load_path filesep 'sta_struct']); sta_struct = sta_struct.sta_struct;
        proc_traces = load([this_load_path filesep 'proc_traces']); proc_traces = proc_traces.proc_traces;
        event_frames = load([this_load_path filesep 'event_frames']); event_frames = event_frames.event_frames;
    else
        %% search for session data in matfile folder
        session_data = load([matfile_dir filesep this_mat_file]);
        % discard session data after session stop time (specified in Tprocc)
        if (~isempty(this_TDTduration)&& ~isnan(this_TDTduration))
            this_tot_frames = round(session_data.sampling_frequency*this_TDTduration);
            session_data.Data_stream_405 = session_data.Data_stream_405(1:this_tot_frames);
            session_data.Data_stream_560 = session_data.Data_stream_560(1:this_tot_frames);
            session_data.Data_stream_465 = session_data.Data_stream_465(1:this_tot_frames);
        end
        
        disp(['loaded single sesison data: ' this_mat_file])
        this_save_path = this_load_path;
        if ~exist(this_save_path, 'dir')
            mkdir(this_save_path)
        end
        %% parameters
        opt = set_opt0(session_data);
        opt.mouse = this_animal;
        opt.date = this_save_date;
        opt.nosepoke_required = 'Nosepoke_L';
        opt.session_file_name = this_mat_file;
        opt.this_save_path = this_save_path;
        opt.sta_pre_sec = 8;
        opt.sta_post_sec = 8;
        opt.baseline_win_start_sec = 8;
        opt.baseline_duration_sec = 4;
        opt.pre_frames = round(opt.sta_pre_sec*opt.frame_rate);
        opt.post_frames = round(opt.sta_post_sec*opt.frame_rate);
        opt.baseline_frames = opt.pre_frames - round(opt.baseline_win_start_sec*opt.frame_rate) +[1:round(opt.frame_rate*opt.baseline_duration_sec)];
        opt.location = this_program;
        opt.plot_sta_fds = {'uv','green'};
        opt.plot_sta_names = {'Isosbestic','Signal'};
        save([this_save_path filesep 'opt0'],'opt')
        disp('running single session analysis')
        [proc_traces,sta_struct,~,event_frames,opt] = run_single_session(session_data,this_save_path,opt,...
            'behav_timestamps',behav_timestamps,...
            'event_fds',event_types,...
            'plot_fds',event_types,'IF_SKIP_BEHAVIOR',1);
        save([this_save_path filesep 'behav_timestamps'],'behav_timestamps')
    end
    %% plot trace with timestamps
    figure('units','normalized','outerposition',[0 0 1, 0.3]); hold on
    if ~isempty(save_struct.SocialOnset)
        arrayfun(@(x,y)plot([x,y],[-3 -3],'color',colors.SocialOnset,'linewidth',1.5),save_struct.SocialOnset,save_struct.SocialOffset)
    end
    
    xlabel('Time (seconds)')
    title([this_animal ' ' this_save_date ' ' this_program])
    ylabel('debleached F, zscored')
    
    % plot photometry trace
    this_trace = downsample(proc_traces.db_green,100);
    plot([1:length(this_trace)]./opt.frame_rate*100,this_trace,'color','black');
    
    % plot session start time (when the intruder was put in)
    plot([1 1].*this_start_time,ylim,'color','r')
    export_fig([this_save_path filesep 'FullTraceWithTimestamps' '.png'])
    saveas(gcf,[this_save_path filesep 'FullTraceWithTimestamps' '.fig'])
    
    %% save 
%     output_mat_name = ['TimeStamps_' this_animal '_' this_save_date '_' this_program];
%     save([timestamp_save_dir filesep output_mat_name],'save_struct')
%     disp(['saved as ' output_mat_name])
    
end

%% plot animal average - run social_sta_summary.m



