%% simple sessions summary
% run 'manual_timestamps_analysis' to process TDT matfilese first;
% this script will load the sta traces and make plots for individual
% sessions and animal average
% activity analysis
excel_file = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\aggression_analysisplan_template.xlsx'; % sessions to process
sheet_name = 'BLA';

session_tag = sheet_name;
Tprocess = readtable(excel_file,'Sheet',sheet_name);
animal_names = unique(Tprocess.Mouse); % OTHERWISE, LIST THE ANIMAL NAMES AND THE SCRIPT WILL IGNORE THE OTHERS

matlab_path = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry'; % CHANGE THIS
data_dir = 'C:\Users\zzhang21\Documents\GitHub\SocialPhotometry\example';
timestamp_dir = [data_dir '\ManualTimestamps']; % put manual timestamp (.csv or .xlsx) files here
fig_save_path = [data_dir '\AggressionFigures']; % figures will be saved here
mat_save_path = [data_dir '\TDTProcData']; % processed data will be saved here
timestamp_save_dir = [timestamp_dir '\ManualTimestamps\Timestamp_files']; %
%% prepare matlab
% add matlab script to path
addpath(genpath(matlab_path));
cd(matlab_path)
% setup fonts
set(0,'defaultAxesFontName','Arial')
set(0,'defaultTextFontName','Arial')
set(0,'defaultfigurecolor',[1 1 1],'DefaultTextColor', [0, 0, 0])
set(0,'DefaultAxesColor',[1 1 1])
set(0,'DefaultAxesXcolor', [0, 0, 0], 'DefaultAxesYcolor', [0, 0, 0], 'DefaultAxesZcolor', [0, 0, 0])
set(0,'defaultAxesTickDir','out');
%% global params
global IF_GET_CORRELATION IF_LOAD_EXIST IF_LOAD_PROCTRACE IF_PLOT_RAW IF_ADD_FIELDS
IF_ADD_FIELDS = false; % load and add field to existing data (temporal use)
IF_LOAD_EXIST = false; % load processed data
IF_GET_CORRELATION = false; % correlogram of red and green
IF_LOAD_PROCTRACE = true;
IF_PLOT_RAW = true;
IF_PLOT_SESSIONS = 0;
color = set_colors;
%% set parameters
% full session trace
plot_trace = 'db_green'; plot_trace_color = color.green;
plot_sta_trace = 'green_procf_sta'; title_tag = [session_tag]; % Sensor_SessionType
plot_rawsta_trace = 'green_procf_sta';
plot_refsta_trace = 'uv_procf_sta'; title_tag = [session_tag]; % Sensor_SessionType

% events if DLC is provided:
plot_events = {'AttackOnset','AttackOffset','AggressionOnset','AggressionOffset','SocialNoAttackOnset','SocialNoAttackOffset'};
plot_event_names = {'Attack onset', 'Attack offset','Aggression onset', 'Aggression offset','Social onset', 'Social offset'};
plot_behav_names = {'Attack','Aggression','SocialNoAttack'};
plot_event_colors = {color.AttackOnset, color.AttackOffset,...
    color.AttackOnset, color.AttackOffset,...
    color.SocialNoAttackOnset, color.SocialNoAttackOffset};

y_label = 'zscored F (debleached)';
plot_frame_rate = 113; %Hz
video_frame_rate = 20;

%% initiate structures
if ~exist(fig_save_path, 'dir')
    mkdir(fig_save_path)
end
session_avg_stas = struct(); session_raw_stas = struct();  session_avg_refstas = struct();
session_episode_mean = struct();
all_stas = struct(); all_raw_stas = struct(); all_ref_stas = struct();
all_peak_times_from_onset = struct(); 
all_peak_times_to_offset = struct(); 
all_fract_time_to_peak = struct();
all_rise_times = struct();
session_animals = struct(); session_count = struct();
for ee = 1:numel(plot_events)
    plot_event = plot_events{ee};
    session_avg_stas.(plot_event) = []; session_raw_stas.(plot_event) = []; session_avg_refstas.(plot_event) = [];
    all_stas.(plot_event) = []; all_raw_stas.(plot_event) = [];  all_ref_stas.(plot_event) = [];
    all_animal_sta_idx.(plot_event) = [];
    all_sta_event_duration.(plot_event) = []; % may differ from all_event_duration
    session_animals.(plot_event) = {}; session_count.(plot_event) = 0;
end
for ee = 1:numel(plot_behav_names)
    plot_event = plot_behav_names{ee};
    all_event_duration.(plot_event) = [];
    session_event_duration.(plot_event) = [];
    all_peak_times_from_onset.(plot_event) = [];
    all_peak_times_to_offset.(plot_event) = [];
    all_fract_time_to_peak.(plot_event) = [];
    all_peak_slope.(plot_event) = [];
    all_peak_activity.(plot_event) = [];
    all_auc.(plot_event) = []; % added for revision
    all_down_slope.(plot_event) = []; % added for revision
    all_animal_dur_idx.(plot_event) = [];
    all_episodes.(plot_event) = [];
end
%% process or load session data
session_counter = 0; session_animal_names = {};

for m = 1:numel(animal_names)
    this_animal = animal_names{m};
    this_rows = find(contains(Tprocess.Mouse,this_animal));
    this_dates = Tprocess.Date(this_rows);
    this_programs =  Tprocess.Program(this_rows);
    for d = 1:numel(this_rows)
        session_counter = session_counter+1;
        session_animal_names{session_counter} = this_animal;
        this_date = datestr(this_dates(d),'yyyymmdd');
        this_program = this_programs{d};
        this_load_path = [mat_save_path filesep this_animal filesep this_date filesep this_program]; % photometry data
        this_timstamp_path = [timestamp_save_dir filesep 'TimeStamps_' this_animal '_' this_date '_' this_program '.mat'];
        
        sta_struct = load([this_load_path filesep 'sta_struct']); sta_struct = sta_struct.sta_struct;
        proc_traces = load([this_load_path filesep 'proc_traces']); proc_traces = proc_traces.proc_traces;
        event_times = load([this_load_path filesep 'behav_timestamps']); event_times = event_times.behav_timestamps;
        plot_event_frames = structfun(@(x)x.*plot_frame_rate,event_times,'un',0);
        opt = load([this_load_path filesep 'opt0']); opt = opt.opt;
        
        % get peak time during each event
        this_trace = proc_traces.(plot_trace);
        trace_event_frames = structfun(@(x)round(x.*opt.frame_rate),event_times,'un',0);
        for ee = 1:numel(plot_behav_names)
            this_behav = plot_behav_names{ee};
            this_onset_frames = trace_event_frames.([this_behav 'Onset']);
            this_offset_frames = trace_event_frames.([this_behav 'Offset']);
            
            % get trial durations for each STA trace
            this_sta_onsets_trials = find(this_onset_frames>opt.pre_frames&this_onset_frames<numel(this_trace)-opt.post_frames); 
            all_sta_event_duration.([this_behav 'Onset']) = [all_sta_event_duration.([this_behav 'Onset']);...
                (this_offset_frames(this_sta_onsets_trials)-this_onset_frames(this_sta_onsets_trials))./opt.frame_rate];
            this_sta_offsets_trials = find(this_offset_frames>opt.pre_frames&this_offset_frames<numel(this_trace)-opt.post_frames); 

            all_sta_event_duration.([this_behav 'Offset']) = [all_sta_event_duration.([this_behav 'Offset']);...
                (this_offset_frames(this_sta_offsets_trials)-this_onset_frames(this_sta_offsets_trials))./opt.frame_rate];
    
            discard_last_events = find(this_onset_frames<opt.pre_frames|this_offset_frames>numel(this_trace)-opt.post_frames);
            this_onset_frames(discard_last_events) = [];
            this_offset_frames(discard_last_events) = [];
            
            % get traces for each event episodes
            all_episodes.(this_behav) = [all_episodes.(this_behav);arrayfun(@(x,y)this_trace(x:y),...
                this_onset_frames,this_offset_frames,'un',0)];
            
            % peak time from onsets
            % add one sec to offset to capture peak
            pad_time = opt.frame_rate;
            this_peak_time_from_onsets = arrayfun(@(x,y)find(smooth(this_trace(x:y),.2*opt.frame_rate)==max(smooth(this_trace(x:y),.2*opt.frame_rate)),1,'first'),...
                this_onset_frames,this_offset_frames);
            this_peak_time_from_offsets = (this_offset_frames-this_onset_frames)-this_peak_time_from_onsets+1;
            this_peak_activity = arrayfun(@(x,y)max(smooth(this_trace(x:y),.2*opt.frame_rate)),...
                this_onset_frames,this_offset_frames)-arrayfun(@(x,y)min(smooth(this_trace(x:x+y-1),.2*opt.frame_rate)),...
                this_onset_frames,this_peak_time_from_onsets);
            this_raw_peak_activity = arrayfun(@(x,y)max(smooth(this_trace(x:y),.2*opt.frame_rate)),...
                this_onset_frames,this_offset_frames);
%             this_peak_downslope_activity = arrayfun(@(x,y)max(smooth(this_trace(x:y),.2*opt.frame_rate)),...
%                 this_onset_frames,this_offset_frames)-arrayfun(@(x,y)min(smooth(this_trace(x-y+1:x),.2*opt.frame_rate)),...
%                 this_offset_frames,this_peak_time_from_offsets);

            this_auc  = arrayfun(@(x,y)trapz(this_trace(x+find(this_trace(x:y)>0))),...
                this_onset_frames,this_offset_frames);
%             this_min_time_from_onsets =  arrayfun(@(x,y)find(smooth(this_trace(x:x+y),.2*opt.frame_rate)==min(smooth(this_trace(x:x+y),.2*opt.frame_rate)),1,'last'),...
%                 this_onset_frames,this_peak_time_from_onsets)./opt.frame_rate;
            this_min_time_from_onsets = 0;
            this_peak_time_from_onsets = this_peak_time_from_onsets./opt.frame_rate;
            this_peak_slope = this_peak_activity./(this_peak_time_from_onsets-this_min_time_from_onsets);
            this_peak_time_to_offsets = -this_peak_time_from_onsets+(this_offset_frames-this_onset_frames)./opt.frame_rate;
%             this_down_slope = this_peak_downslope_activity./this_peak_time_to_offsets;
            all_peak_times_from_onset.(this_behav)=[ all_peak_times_from_onset.(this_behav);this_peak_time_from_onsets(:)];
            all_peak_slope.(this_behav)=[ all_peak_slope.(this_behav);this_peak_slope(:)];
%             all_down_slope.(this_behav)=[ all_down_slope.(this_behav);this_down_slope(:)];
            all_auc.(this_behav)=[ all_auc.(this_behav);this_auc(:)];
            all_peak_activity.(this_behav)=[ all_peak_activity.(this_behav);this_raw_peak_activity(:)];
            all_peak_times_to_offset.(this_behav)=[ all_peak_times_to_offset.(this_behav);this_peak_time_to_offsets(:)];
            all_fract_time_to_peak.(this_behav)= [all_fract_time_to_peak.(this_behav);this_peak_time_from_onsets./((this_offset_frames-this_onset_frames+1)./opt.frame_rate)];
            all_event_duration.(this_behav) = [all_event_duration.(this_behav);...
               (this_offset_frames-this_onset_frames)./opt.frame_rate];
             all_animal_dur_idx.(this_behav) = [all_animal_dur_idx.(this_behav);  ones(length(this_onset_frames),1).*m];

           if numel(this_offset_frames)~=numel(this_peak_slope)
               pause
           end
        end
       

        %% get event duration, in seconds
         for ee = 1:numel(plot_behav_names)
            plot_event = plot_behav_names{ee};
            session_event_duration.(plot_event) = [session_event_duration.(plot_event);...
            mean(event_times.([plot_event 'Offset'])-event_times.([plot_event 'Onset']))];
         
         end
        %% save to animal sta
        for ee = 1:numel(plot_events)
            plot_event = plot_events{ee};
            plot_event_color = plot_event_colors{ee};
            plot_event_name = plot_event_names{ee};
            this_sta_traces = sta_struct.(plot_sta_trace).(plot_event);
            this_refsta_traces = sta_struct.(plot_refsta_trace).(plot_event);

            if isempty(this_sta_traces)
                continue
            else
                session_count.(plot_event) = session_count.(plot_event)+1;
                session_animals.(plot_event){session_count.(plot_event)} = this_animal;
                
                session_avg_stas.(plot_event) = [session_avg_stas.(plot_event);mean(this_sta_traces,1)];
                session_avg_refstas.(plot_event) = [session_avg_stas.(plot_event);mean(this_refsta_traces,1)];
                session_raw_stas.(plot_event) = [session_raw_stas.(plot_event);mean(sta_struct.(plot_rawsta_trace).(plot_event),1)];
                all_stas.(plot_event) = [all_stas.(plot_event); this_sta_traces];
                all_ref_stas.(plot_event) = [all_ref_stas.(plot_event); this_refsta_traces];

                all_animal_sta_idx.(plot_event) = [all_animal_sta_idx.(plot_event); ones(size(this_sta_traces,1),1).*m];
                all_raw_stas.(plot_event) = [all_raw_stas.(plot_event); sta_struct.(plot_rawsta_trace).(plot_event)];
                
            end
            %% plot session trace with behavioral event of interest
            if IF_PLOT_SESSIONS
                figure('name','debleached trace with events','units','normalized','outerposition',[0 0 1 .3]);
                hold on
                ds_factor = round(opt.frame_rate/plot_frame_rate); % downsample to plot_frame_rate Hz
                this_plot_trace = downsample(medfilt1(proc_traces.(plot_trace),10), ds_factor);
                this_plot_x = (1:length(this_plot_trace))./plot_frame_rate;
                plot(this_plot_x, this_plot_trace,'color',plot_trace_color);
                ylabel(y_label); xlabel('Time from session start (sec)');
                plot(xlim,[0 0],'color',[.5 .5 .5]);
                title(['Full session, ' num2str(round(opt.num_frames/opt.frame_rate)), ' sec, display rate ' num2str(plot_frame_rate) ' Hz' ])
                arrayfun(@(x)plot([x x]./opt.frame_rate,ylim,'color',plot_event_color),plot_event_frames.(plot_event));
                this_title = [this_animal,'-',this_date, '_FullTrace_' title_tag];
                suptitle(strrep(this_title,'_',' '))
                export_fig([fig_save_path filesep this_title], '-pdf','-pdf','-painters')
                
                %% Individual trials, raw traces and shaded errorbar
                h=figure('name','sta','units','normalized','outerposition',[0 0 .7 1]);
                ylimit = [-5 5]; y_label = 'debleached F';
                
                subplot(1,3,1); hold on
                imagesc(this_sta_traces);colormap(b2r(ylimit(1),ylimit(2)));
                plot([1 1].*(opt.pre_frames+1),xlim,'color',[.5 .5 .5])
                plot([0 opt.frame_rate]-.5,[0 0],'color','black','linewidth',2);
                text(0,0,{'1 sec'},'color','black','verticalalignment','top');
                colorbar('eastoutside')
                ylim([0 size(this_sta_traces,1)]); xlim([1, size(this_sta_traces,2)])
                ylabel('Trials');
                set(gca,'xtick',[],'xcolor','w')
                xlabel(['Time from ' plot_event_name ' (sec)']);
                axis square
                
                subplot(1,3,2); hold on
                num_trials = size(this_sta_traces,1);
                this_xticks = (-opt.pre_frames:1:0:1:opt.post_frames)./opt.frame_rate;
                arrayfun(@(x)plot(this_xticks,this_sta_traces(x,:)+ylimit(2).*x,'color','black'),1:num_trials);
                plot([0 0],ylim,'color',plot_event_color,'linewidth',2);
                plot(-10.*[1,1],[0,ylimit(2)],'color','black','linewidth',2);
                text(-11,0,{'scale bar';'.5 debleached F'},'color','black','verticalalignment','bottom','horizontalalignment','right');
                axis square
                ylabel('Trials');    set(gca,'ytick',[],'ycolor','w')
                xlabel(['Time from ' plot_event_name ' (sec)']);
                title(['Trial 1-' num2str(num_trials) ', bottom-top'])
                
                subplot(1,3,3); hold on
                shadedErrorBar(this_xticks,nanmean(this_sta_traces,1), nanstd(this_sta_traces,[],1),{'color',plot_trace_color,'linewidth',2},0);
                ylim(ylimit)
                plot([0 0],ylimit,'color',plot_event_color,'linewidth',2);
                plot(xlim,[0 0],'color','black')
                box off; axis square;
                xlabel(['Time from ' plot_event_name ' (sec)']);
                ylabel(y_label)
                axis square
                title('Trial-average')
                this_title = [this_animal,'-',this_date, '_STA_' title_tag];
                suptitle(strrep(this_title,'_',' '))
                export_fig([fig_save_path filesep this_title], '-pdf','-pdf','-painters')
            end
        end
        %     print('-dpdf','-r300', '-painters','myVectorFile2.pdf')%uggly: has large white borders
        %     set(gcf,'renderer','Painters')
        %   exportgraphics(h,[fig_save_path filesep this_title,'.pdf'],'ContentType','vector')
    end
end
%% peak time during events, fraction from start
%% aggression and attack duration and signal peak times
for ee = 1:numel(plot_behav_names)
    plot_event = plot_behav_names{ee};
    % average of all trials
    mouse_rows = all_animal_dur_idx.([plot_event]);
    row_groups = arrayfun(@(x)find(mouse_rows==x),unique(mouse_rows),'un',0);
%     all_fract_time_to_peak.(plot_event) = all_peak_times_from_onset.(plot_event)./(all_peak_times_from_onset.(plot_event)+...
%         all_peak_times_to_offset.(plot_event));
%     
    animal_event_duration.(plot_event) = mean_of_rows(all_event_duration.(plot_event),row_groups);
    animal_peak_times_from_onset.(plot_event) = mean_of_rows(all_peak_times_from_onset.(plot_event),row_groups);
    animal_peak_times_to_offset.(plot_event) = mean_of_rows(all_peak_times_to_offset.(plot_event),row_groups);
    animal_fract_time_to_peak.(plot_event) = mean_of_rows(all_fract_time_to_peak.(plot_event),row_groups);
    
end
%% get animal average from session average
animal_avg_stas = struct();
animal_avg_refstas = struct();
for ee = 1:numel(plot_events)
    plot_event = plot_events{ee};
    % average of session average
    %     [valid_animal_names,~, mouse_rows]= unique(session_animals.(plot_event));
    %     row_groups = arrayfun(@(x)find(mouse_rows==x),unique(mouse_rows),'un',0);
    %     animal_avg_stas.(plot_event) = mean_of_rows(session_avg_stas.(plot_event),row_groups);
    
    % average of all trials
    mouse_rows = all_animal_sta_idx.(plot_event);
    row_groups = arrayfun(@(x)find(mouse_rows==x),unique(mouse_rows),'un',0);
    animal_avg_stas.(plot_event) = mean_of_rows(all_stas.(plot_event),row_groups);
    animal_avg_refstas.(plot_event) = mean_of_rows(all_ref_stas.(plot_event),row_groups);

end


%% save animal average sta
save([mat_save_path filesep 'animal_avg_stas' title_tag],'animal_avg_stas')
save([mat_save_path filesep 'all_stas' title_tag],'all_stas')
save([mat_save_path filesep 'animal_event_duration' title_tag],'animal_event_duration')
save([mat_save_path filesep 'animal_peak_times_from_onset' title_tag],'animal_peak_times_from_onset')
save([mat_save_path filesep 'animal_peak_times_to_offset' title_tag],'animal_peak_times_to_offset')
save([mat_save_path filesep 'animal_fract_time_to_peak' title_tag],'animal_fract_time_to_peak')
save([mat_save_path filesep 'all_fract_time_to_peak' title_tag],'all_fract_time_to_peak')

disp('session avg data saved')
%% sessions mean
pre_avg_sec = 2;
post_avg_sec = 2;
ylimit = [-2,1.5];
xlimit = [-3 3];
bin_frames = {opt.pre_frames+[-pre_avg_sec*opt.frame_rate:-1],...
    opt.pre_frames+[1:post_avg_sec*opt.frame_rate]};
bin_names = {'Pre',...
    'Post' };
session_avg = struct(); % log pre post values
for ee = 1:numel(plot_events)
    plot_event = plot_events{ee};
    plot_event_color = plot_event_colors{ee};
    figure('name','sta','units','normalized','outerposition',[0 0 .7 1]);
    % pooling all trials from all sessions
    num_cols = 3;
    this_sta_traces = all_stas.(plot_event);
    this_refsta_traces = all_ref_stas.(plot_event);

    subplot(2,num_cols,1); hold on
    num_trials = size(this_sta_traces,1);
    this_xticks = (-opt.pre_frames:1:0:1:opt.post_frames)./opt.frame_rate;
    this_frame_range = [(opt.pre_frames+opt.frame_rate*xlimit(1)):(opt.pre_frames+opt.frame_rate*xlimit(2))];
    imagesc(this_sta_traces(:,this_frame_range)); colormap(b2r(-3,3)); colorbar('location','southoutside');
    plot((opt.frame_rate*abs(xlimit(1))+1)*[1,1],[0 num_trials+0.5],'color','black','linewidth',1.5);
    
    % arrayfun(@(x)plot(this_xticks,this_sta_traces(x,:)+x*ylimit(2),'color','black'),1:num_trials);
    % plot([0 0],ylim,'color',plot_event_color,'linewidth',2);
    % plot(-10.*[1,1],[0,ylimit(2)],'color','black','linewidth',2);
    % text(-11,0,{'scale bar';'.5 \DeltaF/F'},'color','black','verticalalignment','bottom','horizontalalignment','right');
    axis square; axis off;
    ylabel('Trials');    set(gca,'ytick',[],'ycolor','w')
    xlabel(['Time from ' plot_event ' (s)']);
    title(['Trial 1-' num2str(num_trials) ', bottom-top'])
    
    subplot(2,num_cols,2); hold on
    shadedErrorBar(this_xticks,nanmean(this_refsta_traces,1), nanstd(this_refsta_traces,[],1)./sqrt(size(this_refsta_traces,1)),{'color',[.5 .5 .5],'linewidth',2},0);
    shadedErrorBar(this_xticks,nanmean(this_sta_traces,1), nanstd(this_sta_traces,[],1)./sqrt(size(this_sta_traces,1)),{'color',plot_trace_color,'linewidth',2},0);
    ylim(ylimit); xlim(xlimit)
    plot([0 0],ylimit,'color',plot_event_color,'linewidth',2);
    plot(xlim,[0 0],'color','black')
    box off; axis square;
    xlabel(['Time from ' plot_event ' (s)']);
    axis square
    title('Trial-average')
    
    subplot(2,num_cols,3); hold on
    this_val = struct();
    for b = 1:numel(bin_names)
        this_val.(bin_names{b}) = nanmean(all_raw_stas.(plot_event)(:,bin_frames{b}),2);
    end
    scatter_cmp_conditions(this_val,[],1,[],'BriefXlabel',1,'connect_scatter',1,...
        'test_type','signrank','plot_stats',1);
    axis square; xtickangle(45); %ylim([-4 4])
    title('Trial-average'); xlabel([]); ylabel('zscored F')
    
    this_sta_traces = animal_avg_stas.(plot_event);
    this_refsta_traces = animal_avg_refstas.(plot_event);

    subplot(2,num_cols,4); hold on
    num_trials = size(this_sta_traces,1);
    this_xticks = (-opt.pre_frames:1:0:1:opt.post_frames)./opt.frame_rate;
    arrayfun(@(x)plot(this_xticks,this_sta_traces(x,:)+x*ylimit(2),'color','black'),1:num_trials);
    plot([0 0],ylim,'color',plot_event_color,'linewidth',2);
    plot(-10.*[1,1],[0,ylimit(2)],'color','black','linewidth',2);
    text(-11,0,{'scale bar'; [ylimit(2) ' Fzscore']},'color','black','verticalalignment','bottom','horizontalalignment','right');
    axis square
    ylabel('Mice');    set(gca,'ytick',[],'ycolor','w')
    xlabel(['Time from ' plot_event ' (s)']);
    title(['Mouse 1-' num2str(num_trials) ', bottom-top'])
    
    subplot(2,num_cols,5); hold on
    shadedErrorBar(this_xticks,nanmean(this_refsta_traces,1), nanstd(this_refsta_traces,[],1)./sqrt(size(this_refsta_traces,1)),{'color',[.5 .5 .5],'linewidth',2},0);
    shadedErrorBar(this_xticks,nanmean(this_sta_traces,1), nanstd(this_sta_traces,[],1)./sqrt(size(this_sta_traces,1)),{'color',plot_trace_color,'linewidth',2},0);
    ylim(ylimit); xlim([-3 3])
    plot([0 0],ylimit,'color',plot_event_color,'linewidth',2);
    plot(xlim,[0 0],'color','black')
    box off; axis square;
    xlabel(['Time from ' plot_event ' (s)']);
    ylabel('Debleached F')
    axis square
    title({'Mouse-average';' '})
    this_title = ['Summary_STA_' plot_event '_' title_tag];
    
    subplot(2,num_cols,6); hold on
    this_val = struct();
    for b = 1:numel(bin_names)
        this_val.(bin_names{b}) = nanmean(this_sta_traces(:,bin_frames{b}),2);
        session_avg.([plot_event '_' bin_names{b}]) = this_val.(bin_names{b});
    end
    scatter_cmp_conditions(this_val,[],1,[],'BriefXlabel',1,'connect_scatter',1,...
        'test_type','signrank','plot_stats',1);
    axis square; xtickangle(45); ylim([-3 3])
    title({'Mouse-average',' '}); xlabel([]); ylabel('zscored F')
     
    suptitle(strrep(this_title,'_',' '))
    exportgraphics(gcf,[fig_save_path filesep this_title '.pdf'],'ContentType','vector')
end

% normalise to pre onset avg.
fds = fields(session_avg);
for f = 1:numel(fds)
    fd = fds{f};
    this_ref_fd = strrep(fd,'Offset','Onset');
    this_ref_fd = strrep(this_ref_fd,'_Post','_Pre');
    session_avg.(['d' fd]) = session_avg.(fd) - session_avg.(this_ref_fd);
end
save([mat_save_path filesep 'session_avg' title_tag],'session_avg')

%% compare pre-onset, post-onset and pre-offset
plot_bin_names = {'Onset_Pre','Onset_Post','Offset_Pre','Offset_Post'};
numplot = numel(plot_behav_names);

figure('name','sta','units','normalized','outerposition',[0 0 numplot*0.2 .8]);
for ee = 1:numel(plot_behav_names)
    plot_behav = plot_behav_names{ee};
    this_val = struct();
    for b = 1:numel(plot_bin_names)
        this_bin_name = plot_bin_names{b};
        this_val_name = [plot_behav,this_bin_name];
        if ~isfield(session_avg,this_val_name)
            continue
        else
            this_val.(this_val_name) = session_avg.(this_val_name);
        end
    end
    subplot(1,numplot,ee)
    scatter_cmp_conditions(this_val,[],1,[],'BriefXlabel',1,'connect_scatter',1,...
        'test_type','signrank','plot_stats',1,'plot_se',1);

    axis square; xtickangle(45); ylim(ylimit)
    ylabel('Fzscored')
    
end
this_title = ['STATS_STA_' title_tag];
suptitle(strrep(this_title,'_',' '))
exportgraphics(gcf,[fig_save_path filesep this_title '.pdf'],'ContentType','vector')

%% compare aggression onset, attack and aggression offset
plot_fd_names = {'AggressionOnset_Pre','AggressionOnset_Post',...
   'AttackOnset_Post',...
   'AggressionOffset_Pre','AggressionOffset_Post'};
% ylimit = [-0.5 1.2];
this_val = struct();
for ee = 1:numel(plot_fd_names)
    this_val_name = plot_fd_names{ee};
    if ~isfield(session_avg,this_val_name)
        continue
    else
        this_val.(this_val_name) = session_avg.(this_val_name);
    end
       
end
figure('name','episode mean','units','normalized','outerposition',[0 0 0.5 .8]);
scatter_cmp_conditions(this_val,[],1,[],'BriefXlabel',1,'connect_scatter',1,...
    'test_type','signrank','plot_stats',1);
[~,p] = structfun(@(x)ttest(x),this_val);
arrayfun(@(x,y)text(x, 1.1*ylimit(2),['p = ' num2str(y,'% 10.3f')],'color',[1 0 0].*(y<0.05),...
    'horizontalalignment','right','verticalalignment','top','Rotation',90),[1:numel(fields(this_val))]',p)
axis square; xtickangle(45); ylim(ylimit)
ylabel('Fzscored')
this_title = ['STATS_AggSTA_' title_tag];
suptitle(strrep(this_title,'_',' '))
exportgraphics(gcf,[fig_save_path filesep this_title '.pdf'],'ContentType','vector')


%% episodes sorted by trial duration
fds = {'Aggression','Attack'};
num_cols = 2;num_rows = numel(fds);plot_count = 1;
figure('name','grouped sta','units','normalized','outerposition',[0 0 num_cols*0.4 1]);
raw_ylimit = [-3,5];
smoothfun = @(x)smooth(x,.2*opt.frame_rate);

for f = 1:numel(fds)
    fd = fds{f};
this_duration = all_event_duration.(fd);
this_episode = all_episodes.(fd);
[~,sorted_idx] = sort(this_duration);

this_mat = zeros(numel(this_duration),max(cellfun(@(x)length(x),this_episode))); % for imagesc plot
for i = 1:numel(this_duration)
%     plot([1:length(this_episode{sorted_idx(i)})]./opt.frame_rate,i+this_episode{sorted_idx(i)},'color','black');
this_mat(i,1:length(this_episode{sorted_idx(i)})) = smoothfun(this_episode{sorted_idx(i)});
end
subplot(num_rows,num_cols,[plot_count,plot_count+num_cols]); hold on
imagesc(this_mat);colormap(b2r(raw_ylimit(1),raw_ylimit(2)))
ylabel(['Sorted ' fd ' episodes'])
xlabel('Frames from trial start')
colorbar('southoutside')
plot_count = plot_count+1;
title(fd)
end

this_title = ['SortedEpisodes' title_tag];
suptitle(strrep(this_title,'_',' '))
exportgraphics(gcf,[fig_save_path filesep this_title '.pdf'],'ContentType','vector')

