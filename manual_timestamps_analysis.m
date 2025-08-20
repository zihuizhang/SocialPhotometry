% 1. read two csv files, one for attack start times and another for attack stop times (e.g. manually labelled from video )
%   - these timestamps are relative to the begining of the video and TDT recording
% 2. read DLC output files (optional), to get social interaction times determined by
% inter-mouse distance
% 3. read TDT mat files to get photometry signal aligned to attack/aggression/social interaction times
% 4. save processed data for each session

% ZZ 2022

% >>>> PLEASE KEEP THE NAMING FORMAT CONSISTENT <<<<<
% - SEPARATE ANIMAL NAME AND DATE BY UNDERSCORES
% - USE THE SAME DATE FORMAT
% - see files names in the example folder
%% set directories
matlab_path = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry'; % CHANGE THIS
excel_file = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\aggression_analysisplan_template.xlsx'; % sessions to process
sheet_name = 'processlist'; % select sheet in excel_file to process

matfile_dir = 'E:\ZoeVideos\TDTMATFiles'; % put TDT converted mat files here. CHANGE THIS
% load_data_dir = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\example'; % CHANGE THIS
load_data_dir = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\example'; % CHANGE THIS
save_data_dir = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\example'; % CHANGE THIS

timestamp_dir = [load_data_dir '\ManualTimestamps']; % put manual timestamp (.csv or .xlsx) files here
dlc_dir = [load_data_dir '\DLCprocdata\']; % put DLC output (.csv) files here
mat_save_path = [save_data_dir '\TDTProcData']; % processed data will be saved here
timestamp_save_dir = [timestamp_dir filesep 'Timestamp_files']; % will save output file to a folder named 'Timestamp_files' under the data directory
addpath(genpath(matlab_path)); cd(matlab_path)
start_fileList = [dir(fullfile(timestamp_dir,'/*Start.csv'));dir(fullfile(timestamp_dir,'/*Start.xlsx'))];
end_fileList = [dir(fullfile(timestamp_dir,'/*End.csv'));dir(fullfile(timestamp_dir,'/*End.xlsx'))];
mat_fileList = dir(fullfile(matfile_dir,'/*.mat'));
if ~exist(timestamp_save_dir,'dir')
    mkdir(timestamp_save_dir);
end
%% set params
video_frame_rate = 20; % Hz; resolution for behav traces
min_attack_dur = 2; % sec used 2 sec for GRAB5HT data
max_attack_isi = 1; % sec
photo_timestamps_name = 'RewardNoShock_timestamps'; % as in the TDT .mat file
event_names = {'Attack'};
colors = set_colors();
Tprocess = readtable(excel_file,'Sheet',sheet_name);
IF_LOAD_PROC = 0;
pixelspermm = 1.7; % for box6.1
dist_thresh = pixelspermm*50; % threshold for defining social interaction
%% Simba excel filelist to search
dlcfileList = dir(fullfile(dlc_dir));

%% global params for photometry
global IF_GET_CORRELATION IF_LOAD_PROCTRACE IF_PLOT_RAW IF_ADD_FIELDS
IF_ADD_FIELDS = false; % load and add field to existing data (temporal use)
IF_GET_CORRELATION = false; % correlogram of red and green
IF_LOAD_PROCTRACE = true;
IF_PLOT_RAW = false;
%% loop through animals
all_start_files = {start_fileList.name};
all_end_files = {end_fileList.name};
all_mat_files = {mat_fileList.name}';
num_sessions = size(Tprocess,1);
for a = 1:num_sessions
    this_animal = Tprocess.Mouse{a};
    this_animal_index_str = strrep(strrep(this_animal,'ZZ',''),'PJA','');
    this_animal_index_pad =  num2str(str2double(this_animal_index_str),'%03.f');
    this_date = Tprocess.Date(a,:);
    this_program = Tprocess.Program{a};
    this_TDTend = Tprocess.TDTend(a,:);

    this_start_time = Tprocess.SessionStartFrame(a,:)/video_frame_rate;
    try
        this_TDTduration = Tprocess.SessionDuration(a,:); % sec
    catch
        this_TDTduration = [];
    end
    photo_duration = Tprocess.PhotoDuration(a,:);
    IF_MANUAL_FRAMES = Tprocess.IF_MANUAL_FRAMES(a,:); % otherwise take as seconds
    
    %% find manual timestamps
    this_manual_date =  datestr(this_date,'mm-dd-yyyy');
    if sum(cellfun(@(x)contains(x,this_animal)&&contains(x,this_manual_date)&&contains(x,this_program),all_start_files))==1
        this_start_file = all_start_files{cellfun(@(x)contains(x,this_animal)&&contains(x,this_manual_date)&&contains(x,this_program),all_start_files)};
    else
        try
            this_manual_date =  datestr(this_date,'yyyymmdd');
            this_start_file = all_start_files{cellfun(@(x)contains(x,this_program)&&contains(x,this_animal)&&contains(x,this_manual_date),all_start_files)};
        catch
            this_start_file = [];
            disp([this_manual_date ' ' this_animal 'not found'])
        end
    end
    if ~isempty(this_start_file)
        this_end_file = strrep(this_start_file,'_Start','_End');
    else
        this_end_file = [];
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
    this_mat_file = all_mat_files{cellfun(@(x)(contains(x,strrep(this_mat_animal,'ZZ','ZZH24-'))||contains(x,strrep(this_mat_animal,'ZZ','ZZH23-'))||contains(x,strrep(this_mat_animal,'ZZ','ZZH22-'))||contains(x,strrep(this_mat_animal,'PJA','PJA24-')))...
        &&contains(x,this_mat_date),all_mat_files)};
    if (~isempty(this_end_file))&&(~isempty(this_mat_file))
        disp(['Processing ' this_mat_file])
    else
        disp(['CHECK INPUT FILES FOR ' this_animal ' ' this_mat_date]);
    end
    this_save_date = datestr(this_date,'yyyymmdd');
    %% get attack and photostim times
    if ~isempty(this_start_file)
        this_attack_start_times = table2array(readtable([timestamp_dir filesep this_start_file]));
        this_attack_end_times = table2array(readtable([timestamp_dir filesep this_end_file]));
    else
        this_attack_start_times = [];this_attack_end_times = [];
    end
    % scale by video recording frame rate if manual stamps were in frames
    if IF_MANUAL_FRAMES
        this_attack_start_times = this_attack_start_times./video_frame_rate;
        this_attack_end_times = this_attack_end_times./video_frame_rate;
    end
    
    % discard events that are too short and merge events close together
    [this_attack_start_times,this_attack_end_times] = filt_event_timestamps(this_attack_start_times,this_attack_end_times,...
        min_attack_dur, max_attack_isi);
    this_attack_durations = this_attack_end_times - this_attack_start_times;

    %% make a data struct same as 'read_plot_excel'
    save_struct = struct();
    save_struct.('AttackOnset') = this_attack_start_times;
    save_struct.('AttackOffset') = this_attack_end_times;
  
    %% load dlc to get aggression episodes 
    % read dlc tables
    input_excel_name = dlcfileList(cellfun(@(x)contains(x,this_animal)&&(contains(x,datestr(Tprocess.Date(a,:),'yyyymmdd'))||contains(x,datestr(Tprocess.Date(a),'mmddyyyy')))...
        &&contains(x,this_program),{dlcfileList(:).name}));
    input_excel_path = [dlc_dir filesep input_excel_name.name];
    mouse_speed = []; mouse_distance = [];
    if (~isempty(input_excel_name))&&exist(input_excel_path,'file')
        Tdata =readtable(input_excel_path); 
        %% get locomotion and distance between mice
        [dlc_vars] = set_dlc_vars('SIMBA_Aggression');
        [mouse_position,~,~,~,...
            mouse_bodyparts,juv_bodyparts] = proc_dlc_data(input_excel_path,'SIMBA_Aggression',[720 540]);
        close
        [distance,interact] = get_intruder_distance(mouse_bodyparts,juv_bodyparts,'dist_thr',dist_thresh,'discard_frames',1:this_start_time*video_frame_rate);
        mouse_speed = arrayfun(@(x)pdist([mouse_position(x,:);mouse_position(x+1,:)],'euclidean'),1:size(mouse_position,1)-1)*video_frame_rate;
        mouse_distance = distance.mouse_juv;
        mouse_speed = mouse_speed./pixelspermm./10;% cm/s
        mouse_distance = mouse_distance./pixelspermm./10;% cm

        % find sniff events with attacks
        this_attack_trace = zeros(size(interact.proximal));
        for i = 1:numel(this_attack_start_times)
            this_attack_trace(round(this_attack_start_times(i)*video_frame_rate):round(this_attack_end_times(i)*video_frame_rate)) = 1;
        end
        this_attack_frames = find(this_attack_trace>0);
        interact.proximal(this_attack_trace>0)=1;
        [event_onset_frames,event_offset_frames] = structfun(@(x)get_onset_frames(find(x>0)'),interact,'UniformOutput',false);
        % discard events that are too short and merge events close together
        [event_onset_frames.proximal,event_offset_frames.proximal] = filt_event_timestamps(event_onset_frames.proximal,event_offset_frames.proximal,...
            min_attack_dur*video_frame_rate, max_attack_isi*video_frame_rate);
        sniff_attack_idx = arrayfun(@(e)(~isempty(intersect(this_attack_frames, event_onset_frames.proximal(e):event_offset_frames.proximal(e)))),1:numel(event_onset_frames.proximal));
        sniff_noattack_idx = ~sniff_attack_idx;
        save_struct.('AggressionOnset') = event_onset_frames.proximal(sniff_attack_idx)./video_frame_rate;
        save_struct.('AggressionOffset') = event_offset_frames.proximal(sniff_attack_idx)./video_frame_rate;
        save_struct.('SocialNoAttackOnset') = event_onset_frames.proximal(sniff_noattack_idx)./video_frame_rate;
        save_struct.('SocialNoAttackOffset') = event_offset_frames.proximal(sniff_noattack_idx)./video_frame_rate;
        save_struct.('SocialOnset') = event_onset_frames.proximal./video_frame_rate;
        save_struct.('SocialOffset') = event_offset_frames.proximal./video_frame_rate;
        
    else
        save_struct.('AggressionOnset')=[]; save_struct.('AggressionOffset') = [];
        save_struct.('SocialNoAttackOnset') = [];  save_struct.('SocialNoAttackOffset') = [];
        save_struct.('SocialOnset') = [];  save_struct.('SocialOffset') = [];
        disp(['reading ' input_excel_path ' error! Skipped'])
    end
    behav_timestamps = save_struct;
    save_struct.speed = mouse_speed;% cm/s
    save_struct.distance = mouse_distance;% cm
        
    %% process photometry data
    this_load_path = [mat_save_path filesep this_animal filesep this_save_date filesep this_program];
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
        opt.nosepoke_required = 'Nosepoke_L'; % dummy
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
            'event_fds',{'AttackOnset','AttackOffset',...
            'AggressionOnset','AggressionOffset','SocialNoAttackOnset','SocialNoAttackOffset','SocialOnset','SocialOffset'},...
            'plot_fds',{'AttackOnset','AttackOffset',...
            'AggressionOnset','AggressionOffset','SocialNoAttackOnset','SocialNoAttackOffset','SocialOnset','SocialOffset'},'IF_SKIP_BEHAVIOR',1);
        save([this_save_path filesep 'behav_timestamps'],'behav_timestamps')
    end
    %% plot trace with timestamps
    figure('units','normalized','outerposition',[0 0 1, 0.3]); hold on
    if ~isempty(save_struct.AttackOnset)
        arrayfun(@(x)plot(x.*[1,1],[-5 5],'color',colors.Attack,'linewidth',1.5),save_struct.AttackOnset)
        arrayfun(@(x)plot(x.*[1,1],[-5 5],'color',colors.Attack,'linestyle',':','linewidth',1.5),save_struct.AttackOffset)
    end
    
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
%     plot([1:length(mouse_distance)]./video_frame_rate,mouse_distance,'color',[.5 .5 .5]);
    export_fig([this_save_path filesep 'FullTraceWithTimestamps' '.png'])
   saveas(gcf,[this_save_path filesep 'FullTraceWithTimestamps' '.fig'])
 

end




