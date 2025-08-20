function [proc_traces,sta_struct,behav_struct,event_frames,opt] = run_single_session(session_data,this_save_path,opt,varargin)
behav_timestamps = [];
plot_fds = {'RewardPortEntry','RewardCue','Nosepoke_R','Nosepoke_L','FirstNosepokeAfterReward'};
event_fds = {'RewardCue','RewardPortEntry','FirstNosepokeAfterReward'};
IF_SKIP_BEHAVIOR = 0; % set to one for imaging only sessions
IF_BEHAVIOR_ONLY = 0;
IF_SKIP_DEBLEACHING = 0;
IF_SKIP_STA = 0;
proc_traces = [];
sta_struct = [];
behav_struct = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'behav_timestamps')
        behav_timestamps = varargin{v+1};
    elseif strcmpi(varargin{v},'plot_fds')
        plot_fds = varargin{v+1};
    elseif strcmpi(varargin{v},'event_fds')
        event_fds = varargin{v+1};
    elseif strcmpi(varargin{v},'IF_SKIP_BEHAVIOR')
        IF_SKIP_BEHAVIOR = varargin{v+1};
    elseif strcmpi(varargin{v},'IF_BEHAVIOR_ONLY')
        IF_BEHAVIOR_ONLY = varargin{v+1};
    elseif strcmpi(varargin{v},'IF_SKIP_DEBLEACHING')
        IF_SKIP_DEBLEACHING = varargin{v+1};
    elseif strcmpi(varargin{v},'IF_SKIP_STA')
        IF_SKIP_STA = varargin{v+1};
    end
end

this_mouse = opt.mouse;
this_date = opt.date;
this_location = opt.location;
session_type = opt.session_type;
event_detection_channel = getOr(opt,'event_detection_channel','green'); % detect events in green channel by default
color = set_colors();
num_frames = opt.num_frames;
global IF_GET_CORRELATION IF_LOAD_PROCTRACE IF_PLOT_RAW IF_GET_SPONTANEOUS

%% get time stamps from TDT file
% add first post-reward magazine entry timestamps as 'RewardPortEntry'
% add 'RewardCue', same as 'Reward'.
event_frames = struct();

try
    [this_stamps,error_stamps] = arrayfun(@(x)findOrnan( session_data.Magazine_timestamps,x),session_data.Reward_timestamps);
    session_data.RewardPortEntry_timestamps =session_data.Magazine_timestamps( this_stamps(~error_stamps)); % rewarded magazine entry
    session_data.RewardCue_timestamps = session_data.Reward_timestamps;
    session_data.StimsTimesOnset_timestamps = session_data.RewardNoShock_timestamps;
    
    % first nosepoke after rewarded magazine entry ('back-to-work')
    [this_stamps,error_stamps] = arrayfun(@(x)findOrnan( session_data.([opt.nosepoke_required '_timestamps']),x),session_data.RewardPortEntry_timestamps);
    session_data.FirstNosepokeAfterReward_timestamps =session_data.([opt.nosepoke_required '_timestamps'])( this_stamps(~error_stamps)); % rewarded magazine entry
    color.FirstNosepokeAfterReward = color.([opt.nosepoke_required]).*0.7;
    timestamp_names = {'RewardPortEntry','Nosepoke_R','Nosepoke_L','Magazine','RewardCue','FirstNosepokeAfterReward','RewardNoShock'};
    for i = 1:numel(timestamp_names)
        event_frames.(timestamp_names{i}) = round(opt.frame_rate.*session_data.([timestamp_names{i} '_timestamps']));
    end
    % event_frames.('RewardCue') = round(opt.frame_rate.*session_data.(['RewardCue' '_timestamps']));
end

%% additional timestamps (e.g. time of attack etc)
if ~isempty(behav_timestamps)
    timestamp_names = fields(behav_timestamps);
    for i = 1:numel(timestamp_names)
        event_frames.(timestamp_names{i}) = round(opt.frame_rate.*behav_timestamps.(timestamp_names{i}));
        session_data.([timestamp_names{i} '_timestamps']) = behav_timestamps.(timestamp_names{i});
    end
end
save([this_save_path filesep 'event_frames'],'event_frames')

%% get raw photometry traces
if ~ IF_BEHAVIOR_ONLY
    % downsample after average filter
    func = @(x)downsample(movmean(x,opt.ds_factor),opt.ds_factor);
    func =  @(x)resample(double(x),opt.frame_rate,round(session_data.sampling_frequency));
    this_red = func(session_data.(opt.red_channel_name));
    this_green = func(session_data.(opt.green_channel_name));
    this_uv = func(session_data.(opt.uv_channel_name));
    
    %% show raw data
    if IF_PLOT_RAW
        figure('name','raw_data','units','normalized','outerposition',[0 0 1 .5]); hold on;
        ylimit = [min([this_red(opt.discard_first_frames:end) this_green(opt.discard_first_frames :end)]),...
            max([this_red(opt.discard_first_frames:end) this_green(opt.discard_first_frames :end)])];
        clear h; h = [];
        h(1) = plot([1:length(this_red)]./opt.frame_rate,this_red,'color',color.red,'displayname','Red');
        h(2) = plot([1:length(this_green)]./opt.frame_rate,this_green,'color',color.green,'displayname','Green');
        h(3) = plot([1:length(this_uv)]./opt.frame_rate,this_uv,'color',color.uv,'displayname','UV');
        try
        for i = 1:numel(timestamp_names)
            this_timestamp = timestamp_names{i};
            this_frames = event_frames.(timestamp_names{i});
            if ~isempty(this_frames)
                for f = 1:numel(this_frames)
                    if f == 1
                        h(i+3) = plot(this_frames(f).*[1,1]./opt.frame_rate,ylimit,'color',color.(this_timestamp),'displayname',this_timestamp);
                    else
                        plot(this_frames(f).*[1,1]./opt.frame_rate,ylimit,'color',color.(this_timestamp));
                    end
                end
            else % plot a dummy
                h(i+3) = plot([0 0],ylimit,'color',color.(this_timestamp),'displayname',this_timestamp);
            end
        end
        end
        ylabel('Raw F'); xlabel('Seconds')
        legend(h(:));
        this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_RAW'];
        suptitle(strrep(this_title,'_',' '))
        export_fig([this_save_path filesep this_title])
        export_fig([this_save_path filesep this_title '.png'])
    end
end
%% --- behavior summary ---
if ~ IF_SKIP_BEHAVIOR
%% make event traces
timestamp_names = fields(event_frames);
for i = 1:numel(timestamp_names)
    this_trace =  zeros(1,opt.behav_num_frames);
    this_frames = event_frames.(timestamp_names{i});
    if ~isempty(this_frames)
        this_trace(this_frames) = 1;
    end
    event_traces.(timestamp_names{i}) = this_trace;
end

%% plot event traces
if IF_PLOT_RAW
    figure('name','session behav','units','normalized','outerposition',[0 0 1 .5]); hold on;
    %     plot_fds = {opt.nosepoke_required,'RewardCue','Magazine','FirstNosepokeAfterReward'};
    h = [];
    try
    for i = 1:numel(plot_fds)
        fd = plot_fds{i};
        h(i) = plot([1:length(event_traces.(fd))]./opt.behav_frame_rate,event_traces.(fd),'color',color.(fd),'displayname',fd);
        plot([1:length(event_traces.(fd))]./opt.behav_frame_rate,i+event_traces.(fd),'color',color.(fd))
        plot(xlim,[i i],'color','w','linewidth',3)
    end
    set(gca,'YTick',[],'YColor','w')
    xlabel('Time from session start (sec)')
    legend(h(1:numel(plot_fds)))
    this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_BehavTrace'];
    suptitle(strrep(this_title,'_',' '))
    export_fig([this_save_path filesep this_title '.png'])
    end
end
end
%% make and plot behav around rewards
if ~IF_SKIP_BEHAVIOR
    behav_struct = struct();
    behav_struct.opt = opt;
    plot_behav_fds = {'Magazine','Nosepoke_L','Nosepoke_R'};
    % event_fds = {'RewardCue','RewardPortEntry','FirstNosepokeAfterReward'};
    for ee = 1:numel(event_fds)
        this_event = event_fds{ee};
        if ~isempty(event_frames.(this_event) )
                figure('name','session behav sta','units','normalized','outerposition',[0 0 .8 1]);
                num_plot_rows = 1; num_plot_cols = numel(plot_behav_fds);
            for i = 1:numel(plot_behav_fds)
                fd = plot_behav_fds{i};
                        subplot(num_plot_rows,num_plot_cols,i); hold on
                this_event_frames = find(event_traces.(this_event)==1);
                this_event_frames(this_event_frames>opt.behav_num_frames-opt.behav_post_frames-1)=[];
                this_event_frames(this_event_frames<opt.behav_pre_frames)=[];
                
                [this_sta_traces] = make_sta_traces(event_traces.(fd), this_event_frames,...
                    opt.behav_pre_frames, opt.behav_post_frames);
                        [yy,xx]=find(this_sta_traces>0);
                        scatter([xx-opt.behav_pre_frames]./opt.behav_frame_rate,yy,...
                            'markerfacecolor',color.(fd),'markeredgecolor','none','marker','square')
                        plot([0 0],ylim,'color',color.RewardCue);
                        xlim([-3 3]); %ylim([0 tot_num_reward+1])
                        xlabel(['Time from ' this_event ' (sec)']); ylabel('Trials'); title(fd)
                        axis square
                
                % get first event time after reward
                this_first_post = arrayfun(@(x)(find(this_sta_traces(x,:))-opt.behav_pre_frames),1:size(this_sta_traces,1),'un',false);
                this_first_post = cellfun(@(x)min(x(x>0)),this_first_post,'un',false);
                this_first_post(cellfun('isempty',this_first_post)) = {nan};
                this_first_post = cell2mat(this_first_post)./opt.behav_frame_rate;
                
                % get last event before reward
                this_last_pre = arrayfun(@(x)(find(this_sta_traces(x,:))-opt.behav_pre_frames),1:size(this_sta_traces,1),'un',false);
                this_last_pre = cellfun(@(x)max(x(x<0)),this_last_pre,'un',false);
                this_last_pre(cellfun('isempty',this_last_pre)) = {nan};
                this_last_pre = cell2mat(this_last_pre)./opt.behav_frame_rate;
                
                % save to struct
                behav_struct.(this_event).([fd '_trace']) = this_sta_traces;
                behav_struct.(this_event).(['first_' fd]) = this_first_post;
                behav_struct.(this_event).(['last_' fd]) = this_last_pre;
                
            end
            
                this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_BehavAround' this_event];
                suptitle(strrep(this_title,'_',' '))
                export_fig([this_save_path filesep this_title '.png'])
        end
    end
    %% save behavior struct
    save([this_save_path filesep 'behav_struct'],'behav_struct')
end
if IF_BEHAVIOR_ONLY
    return
end

%% ---  pre-process raw traces ---
if IF_LOAD_PROCTRACE
    try
        proc_traces = load([this_save_path, filesep 'proc_traces']);proc_traces = proc_traces.proc_traces;
    catch
        disp('Load ProcTrace error! Processing raw trace...')
        IF_LOAD_PROCTRACE = false;
    end
end
if ~ IF_LOAD_PROCTRACE
    if ~IF_SKIP_DEBLEACHING
        %% subtract uv without debleaching
        % smooth uv channel before subtraction
%         [this_green] = get_uv_subtracted_trace(this_green,this_uv);
%         [this_red] = get_uv_subtracted_trace(this_red,this_uv);        
        
        %% de-bleaching, Barker 2020
        disp('debleaching...')
        [this_db_green] = debleach(this_green,opt);
        [this_db_red] =  debleach(this_red,opt);
        [this_db_uv] = debleach(this_uv,opt);
        disp('Done.')
%     else

        % plot_raw_vs_proc(zscore(this_red),this_db_red-this_smoothed_uv,color.red,'Processed',opt)
        % this_title = ['SubUVafterDebleach ' this_mouse '_' this_date '_' session_type '_' this_location '_Red'];
        % suptitle(strrep(this_title,'_',' '));
        % export_fig([this_save_path filesep this_title '.png'])
        
    end
    % show detrend trace compared to raw
%     plot_raw_vs_proc(zscore(this_red),this_db_red,color.red,'DeBleached',opt)
%     this_title = ['DebleachCheck ' this_mouse '_' this_date '_' session_type '_' this_location '_Red'];
%     suptitle(strrep(this_title,'_',' '));export_fig([this_save_path filesep this_title '.png'])
    plot_raw_vs_proc(zscore(this_green),this_db_green,color.green,'DeBleached',opt)
    this_title = ['DebleachCheck ' this_mouse '_' this_date '_' session_type '_' this_location '_Green'];
    suptitle(strrep(this_title,'_',' '));export_fig([this_save_path filesep this_title '.png'])
    
    
    %% save processed traces
    proc_traces = struct();
    proc_traces.opt = opt;
    proc_traces.raw_red = this_red;
    proc_traces.raw_green = this_green;
    proc_traces.raw_uv = this_uv;
    proc_traces.db_red = this_db_red;
    proc_traces.db_green = this_db_green;
    proc_traces.db_uv = this_db_uv;
    save([this_save_path filesep 'proc_traces'],'proc_traces')
end
%% --- get event triggered average signal ---
if ~IF_SKIP_STA
sta_struct = struct();
sta_struct.opt = opt;
channels = {'red','green','uv'};
timestamp_names = fields(event_frames);
for c = 1:numel(channels)
    this_channel = channels{c};
    % RAW dF/F
    this_trace = func(session_data.(opt.([this_channel '_channel_name'])));
    for i = 1:numel(timestamp_names)
        this_timestamp = timestamp_names{i};
        this_frames = event_frames.(timestamp_names{i});
        [sta_traces,~,dff_sta_traces] = make_sta_traces(double(this_trace), this_frames,...
            opt.pre_frames, opt.post_frames, 'baseline_frames',opt.baseline_frames);
        %          sta_struct.([this_channel '_sta']).(this_timestamp) = sta_traces;
        sta_struct.([this_channel '_dff_sta']).(this_timestamp) = dff_sta_traces;
        sta_struct.([this_channel '_f_sta']).(this_timestamp) = sta_traces;
        
    end
    
    % zscored, processed F
    this_trace = proc_traces.(['db_' this_channel ]);
    for i = 1:numel(timestamp_names)
        this_timestamp = timestamp_names{i};
        this_frames = event_frames.(timestamp_names{i});
        [sta_traces,df_sta_traces] = make_sta_traces(double(this_trace), this_frames,...
            opt.pre_frames, opt.post_frames, 'baseline_frames',opt.baseline_frames);
        %         sta_struct.([this_channel '_proc_sta']).(this_timestamp) = sta_traces;
        sta_struct.([this_channel '_dproc_sta']).(this_timestamp) = df_sta_traces;
        sta_struct.([this_channel '_procf_sta']).(this_timestamp) = sta_traces;
    end
    
end

%% temp test
% ~ 0.2 delay from the stimulation time stamps to the actual photostim
% start?
% this_frames = round(session_data.([timestamp_names{i} '_timestamps']).*session_data.sampling_frequency);
% this_trace = func(session_data.(opt.(['green' '_channel_name'])));
% pre_frames = round(5*session_data.sampling_frequency);
% post_frames = round(10*session_data.sampling_frequency);
% baseline_frames = 1:pre_frames(end);
%  [temp_sta_traces,~,temp_dff_sta_traces] = make_sta_traces(double(this_trace), this_frames,...
%             pre_frames, post_frames, 'baseline_frames',baseline_frames);
%% save sta struct
save([this_save_path filesep 'sta_struct'],'sta_struct','-v7.3')

%% plot sta traces
% sta from raw F
figure('name','raw F STA','units','normalized','outerposition',[0 0 1 1])
plot_channels = opt.plot_sta_fds;
plot_trace_fd = '_dff_sta'; y_label = '\Delta F/F';
ylimits = struct();
ylimits.green = [-.1 .1]; ylimits.red = [-0.1 0.1]; ylimits.uv = [-0.1 0.1];
plot_sta_traces(sta_struct,plot_trace_fd,plot_channels,plot_fds,ylimits,y_label,color,opt)
this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_STA(fromRaw)'];
suptitle(strrep(this_title,'_',' '))
export_fig([this_save_path filesep this_title '.png'])
% export_fig([this_save_path filesep this_title '.pdf'],'-pdf','-pdf','-painters')

% sta from processed trace (debleach and uv subtraction)
figure('name','proc STA','units','normalized','outerposition',[0 0 1 1])
plot_trace_fd = '_dproc_sta'; y_label = '\Delta zscore';
% _procf_sta
plot_trace_fd = '_procf_sta'; y_label = 'Debleached zscore';

ylimits = struct();
ylimits.green = [-5 5]; ylimits.red = [-5 5]; ylimits.uv = [-5 5];
plot_sta_traces(sta_struct,plot_trace_fd,plot_channels,plot_fds,ylimits,y_label,color,opt)
this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_STA(fromProc)'];
suptitle(strrep(this_title,'_',' '))
export_fig([this_save_path filesep this_title '.png'])
% export_fig([this_save_path filesep this_title '.pdf'],'-pdf','-pdf','-painters')
saveas(gcf,[this_save_path filesep this_title '.fig'])

%% Sorted STA: DA amplitude vs. reward port entry time
% opt.sta_amp_win = [0,1]; % time window to average after event time
% avg_frames = [opt.sta_amp_win(1)*opt.frame_rate : opt.sta_amp_win(end)*opt.frame_rate]+opt.pre_frames;
%
% figure('name','sorted sta','units','normalized','outerposition',[0 0 .5 1]); hold on;
% num_plot_rows = 2; num_plot_cols = numel(plot_channels);
%
% sta_fds = opt.plot_sta_fds;
% sta_names = opt.plot_sta_names;
%
% sta_type = '_dproc_sta';
% sta_unit = '\Delta zscore';
% event_fd = 'RewardCue';
%
% for s = 1:numel(sta_fds)
%     sta_fd = [sta_fds{s} sta_type];
%     sta_traces =  sta_struct.(sta_fd).(event_fd);
%     sta_traces = filt_sta_traces(sta_traces,opt.filt_param,'METHOD','lowess');
%     first_magazine = behav_struct.(event_fd).('first_Magazine');
%     num_trials = min( size(sta_traces,1),numel(first_magazine));
%     sta_traces = sta_traces(1:num_trials,:); first_magazine = first_magazine(1:num_trials);
%     sta_amps = max( sta_traces(:,avg_frames),[],2);
%     [sorted_first_magazine, sorted_trial_idx] = sort(first_magazine);
%     ax = subplot(num_plot_rows,num_plot_cols,[1:num_plot_cols:(num_plot_rows-1)*num_plot_cols]+(s-1));
%     hold on
%     imagesc(sta_traces(sorted_trial_idx,:)); colormap(ax,b2r(-4,4))
%     scatter(sorted_first_magazine.*opt.frame_rate+opt.pre_frames,1:num_trials,...
%         'square','markeredgecolor',color.('RewardPortEntry'),'markerfacecolor','none','markeredgealpha',.5)
%     set(gca,'xcolor','w')
%     plot([0 opt.frame_rate]-1,[0 0],'color','black','linewidth',2);
%     text(0,0,{'1 sec'},'color','black','verticalalignment','top');
%     ylim([0 0.5+num_trials]); ylabel('Sorted trials');
%     xlim([0 size(sta_traces,2)]);
%     cb = colorbar('southoutside','box','off'); set(get(cb,'label'),'string',sta_unit)
%
%     colorbar('southoutside')
%     plot([1 1].*opt.pre_frames,ylim,'color',color.(event_fd),'linewidth',2,'linestyle',':');
%     title({sta_names{s};['aligned to ' event_fd]})
%
%     subplot(num_plot_rows,num_plot_cols,(num_plot_rows-1)*num_plot_cols+s);
%     hold on
%     yy = sta_amps; xx = first_magazine;
%     y_label = sta_names{s}; x_label = 'Reward port entry (sec)';
%     ylimit = [-0.1 .1]; xlimit = [0 3];
%     scatter(xx,yy,10,'markerfacealpha',.9,'markeredgecolor','none','markerfacecolor',[0 0 0],...
%         'markerfacealpha',1),box off; axis square
%     try
%         [~,~,~,~,~,~,p,r] = plot_fit_linear( xx,yy,xlimit,[0 0 0],'linestyle','--');
%         text(1,1,['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')],'units','normalized', ...
%             'horizontalalignment','right','verticalalignment','top', 'color',[0 0 0])
%     end
%     xlabel(x_label);ylabel(y_label);
% end
% this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_SortedSTA'];
% suptitle(strrep(this_title,'_',' '))
% export_fig([this_save_path filesep this_title '.png'])

%% Session summary plot
for i = 1:numel(event_fds)
    fd = event_fds{i};
    if ~isempty(event_frames.([fd ]))
        figure('name','individual trials','units','normalized','outerposition',[0 0 1 1])
        
        plot_trace_fd = '_procf_sta'; % zscored debleached F without baselining
        num_plot_rows = 3; num_plot_cols = numel(plot_channels)+1;
        
        % traces and events
        subplot(num_plot_rows,num_plot_cols,1:num_plot_cols:num_plot_rows*num_plot_cols); hold on;
        plot_step = 4;
        this_xticks = (-opt.pre_frames:1:0:1:opt.post_frames)./opt.frame_rate;
        channels = {opt.plot_sta_fds};
        trial_count = 0;
        for cc = 1:numel(channels)
            this_channels = channels{cc};
            % plot traces from two channels
            for c = 1:numel(this_channels)
                this_channel = this_channels{c};
                this_filt = filt_sta_traces(sta_struct.([this_channel plot_trace_fd]).(fd),opt.filt_param,'METHOD','lowess');
                arrayfun(@(x)plot(this_xticks,this_filt(x,:)+plot_step*(trial_count+x-1),'color',color.(this_channel)),1:size(this_filt,1))
            end
            % mark reward port entry, and poke times
            ffds = plot_fds;
            for ff = 1:numel(ffds)
                ffd = ffds{ff};
                this_stamps = session_data.([ffd '_timestamps']);
                if ~isempty(this_stamps)
                    for x = 1:size(this_filt,1)
                        this_ffd_stamps = this_stamps((this_stamps<session_data.([fd '_timestamps'])(x)+opt.sta_post_sec(end))&...
                            (this_stamps>session_data.([fd '_timestamps'])(x)-opt.sta_pre_sec))-session_data.([fd '_timestamps'])(x);
                        if ~isempty(this_ffd_stamps)
                            scatter(this_ffd_stamps,plot_step*(trial_count+x-1).*ones(size(this_ffd_stamps)),'square','markerfacecolor',color.(ffd),'markeredgecolor','none','markerfacealpha',.7)
                        end
                    end
                end
            end
            trial_count = trial_count+size(this_filt,1);
        end
        xlim([-opt.sta_pre_sec opt.sta_post_sec])
        ylim([-plot_step plot_step*(trial_count)])
        set(gca,'ycolor','w'); plot([0 0],ylim,'color',color.(fd))
        xlabel(['Time from ' fd ' (sec)']);
        title({fd; [num2str(trial_count),' trials']})
        
        % imagesc
        channels = opt.plot_sta_fds;
        plot_channel_names =  opt.plot_sta_names;
        ylimit = [-plot_step,plot_step];
        for cc = 1:numel(channels)
            subplot(num_plot_rows,num_plot_cols,cc+[1:num_plot_cols:num_plot_rows*num_plot_cols]); hold on;
            this_channels = {channels{cc}};
            this_traces = [];
            
            for c = 1:numel(this_channels)
                this_channel = this_channels{c};
                this_filt = filt_sta_traces(sta_struct.([this_channel plot_trace_fd]).(fd),opt.filt_param,'METHOD','lowess');
                arrayfun(@(x)plot(this_xticks,this_filt(x,:)+plot_step*(trial_count+x-1),'color',color.(this_channel)),1:size(this_filt,1))
                this_traces = [this_traces;this_filt];
                % add a dummy trial to mark trial type separation
                this_traces = [this_traces;zeros(1,size(this_traces,2))];
            end
            imagesc(this_traces);colormap(b2r(ylimit(1),ylimit(2)));
            plot([1 1].*(opt.pre_frames+1),xlim,'color',[.5 .5 .5])
            plot([0 opt.frame_rate]-.5,[0 0],'color','black','linewidth',2);
            text(0,0,{'1 sec'},'color','black','verticalalignment','top');
            colorbar('southoutside')
            ylim([0 size(this_traces,1)]); xlim([1, size(this_traces,2)])
            ylabel('Trials'); title(plot_channel_names{cc})
            set(gca,'xtick',[],'xcolor','w')
            xlabel({'';['Time from ' fd ]});
        end
        
        this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_IndividualTrials_' fd];
        suptitle(strrep(this_title,'_',' '))
        export_fig([this_save_path filesep this_title '.png'])
    end
end
end

%% --- correlation between two channels ---
if IF_GET_CORRELATION
    corr_struct = struct();
    corr_struct.opt = opt;
    trace_fd = '_dproc_sta'; y_label = '\Delta zscore';
    event_fds = {'RewardPortEntry','RewardCue'};
    % filt_param = round(0.03*opt.frame_rate); % window size for midandmean filter
    this_xticks = [opt.sta_corr_range(1):opt.sta_corr_range(end)]./opt.frame_rate-opt.sta_pre_sec;
    for f = 1:numel(event_fds)
        fd = event_fds{f};
        % smooth statraces before taking correlation
        this_red_filt = filt_sta_traces(sta_struct.(['red' trace_fd]).(fd),opt.filt_param,'METHOD','lowess');
        this_green_filt = filt_sta_traces(sta_struct.(['green' trace_fd]).(fd),opt.filt_param,'METHOD','lowess');
        % confine frame rang for correlation
        this_red_filt = this_red_filt(:,opt.sta_corr_range(1):opt.sta_corr_range(end));
        this_green_filt = this_green_filt(:,opt.sta_corr_range(1):opt.sta_corr_range(end));
        % correlation between red and green in individual trials
        [shifts,correlation,cshifts,xcenter] = get_row_corr(this_red_filt,this_green_filt);
        corr_lag_time = shifts./opt.frame_rate; % negative means green precedes red
        % peak time difference
        peak_lag_time = arrayfun(@(x)(find(this_red_filt(x,:)==max(this_red_filt(x,:)))-...
            find(this_green_filt(x,:)==max(this_green_filt(x,:)))),1:size(this_red_filt,1),'un',false);
        peak_lag_time(cellfun(@(x)isempty(x),peak_lag_time)) = {nan};
        peak_lag_time = cell2mat(peak_lag_time)./opt.frame_rate;
        corr_struct.(fd).corr_lag_time = corr_lag_time;
        corr_struct.(fd).cshifts = cshifts;
        corr_struct.(fd).xcenter = xcenter;
        corr_struct.(fd).correlation = correlation;
        corr_struct.(fd).peak_lag_time = peak_lag_time;
        
        
        %% plot
        figure('name','red green correlation trialwise','units','normalized','outerposition',[0 0 1 1])
        num_plot_rows = 3; num_plot_cols = 3;
        subplot(num_plot_rows,num_plot_cols,1:num_plot_cols:num_plot_rows*num_plot_cols); hold on;
        plot_step = 4;
        
        % plot traces from two channels
        arrayfun(@(x)plot(this_xticks,this_red_filt(x,:)+plot_step*(x-1),'color',color.red),1:size(this_red_filt,1))
        arrayfun(@(x)plot(this_xticks,this_green_filt(x,:)+plot_step*(x-1),'color',color.green),1:size(this_green_filt,1))
        try
            arrayfun(@(x)scatter(this_xticks(this_red_filt(x,:)==max(this_red_filt(x,:))),max(this_red_filt(x,:))+plot_step*(x-1),'*','markeredgecolor',color.red,'markerfacecolor','none'),1:size(this_green_filt,1))
            arrayfun(@(x)scatter(this_xticks(this_green_filt(x,:)==max(this_green_filt(x,:))),max(this_green_filt(x,:))+plot_step*(x-1),'*','markeredgecolor',color.green,'markerfacecolor','none'),1:size(this_green_filt,1))
        end
        % mark reward port entry, and poke times
        ffds = {'RewardCue','Magazine',opt.nosepoke_required};
        for ff = 1:numel(ffds)
            ffd = ffds{ff};
            this_stamps = session_data.([ffd '_timestamps']);
            for x = 1:size(this_green_filt,1)
                this_ffd_stamps = this_stamps((this_stamps<session_data.([fd '_timestamps'])(x)+opt.sta_corr_range_sec(end))&...
                    (this_stamps>session_data.([fd '_timestamps'])(x)+opt.sta_corr_range_sec(1)))-session_data.([fd '_timestamps'])(x);
                if ~isempty(this_ffd_stamps)
                    scatter(this_ffd_stamps,plot_step*(x-1).*ones(size(this_ffd_stamps)),'square','markerfacecolor',color.(ffd),'markeredgecolor','none','markerfacealpha',.7)
                end
            end
        end
        ylim([-plot_step plot_step*(x)])
        set(gca,'ycolor','w'); plot([0 0],ylim,'color',color.(fd))
        xlabel('Time from event (sec)')
        title({fd; [num2str(x),' trials']})
        
        subplot(num_plot_rows,num_plot_cols,[1:num_plot_cols:num_plot_rows*num_plot_cols]+1);hold on
        imagesc(correlation); colormap(b2r(-1,1));
        set(gca, 'YDir','normal')
        xlim([1,size(correlation,2)])
        scatter(cshifts,1:size(correlation,1),'*','markeredgecolor','black','markerfacecolor','none')
        % mark reward port entry, and poke times
        ffds = {'RewardCue','Magazine',opt.nosepoke_required};
        for ff = 1:numel(ffds)
            ffd = ffds{ff};
            this_stamps = session_data.([ffd '_timestamps']);
            for x = 1:size(this_green_filt,1)
                this_ffd_stamps = (this_stamps((this_stamps<session_data.([fd '_timestamps'])(x)+opt.sta_corr_range_sec(end))&...
                    (this_stamps>session_data.([fd '_timestamps'])(x)+opt.sta_corr_range_sec(1)))-session_data.([fd '_timestamps'])(x)).*opt.frame_rate+xcenter;
                if ~isempty(this_ffd_stamps)
                    scatter(this_ffd_stamps,x.*ones(size(this_ffd_stamps)),'square','markerfacecolor',color.(ffd),'markeredgecolor','none','markerfacealpha',.7)
                end
            end
        end
        ylim([0 x+1]); ylabel('Trials')
        plot(xcenter.*[1,1],ylim,'color',[.5 .5 .5])
        xlabel(['Frames around ' fd])
        title(['Correlogram Red&Green'])
        
        subplot(num_plot_rows,num_plot_cols,3); hold on
        shadedErrorBar(this_xticks,nanmean(this_red_filt,1), nanstd(this_red_filt,[],1),{'color',color.red,'linewidth',2},0.1);
        shadedErrorBar(this_xticks,nanmean(this_green_filt,1), nanstd(this_green_filt,[],1),{'color',color.green,'linewidth',2},0.1);
        plot([0 0],ylim,'color',color.(fd),'linewidth',2,'linestyle',':');
        box off; axis square;
        xlabel('Time from event (sec)'); ylabel(y_label)
        
        % average trace peak and correlation lag
        subplot(num_plot_rows,num_plot_cols,3+num_plot_cols); hold on
        histogram(corr_lag_time,'binwidth',.5,'facecolor','black','edgecolor','none')
        plot(median(corr_lag_time).*[1,1],ylim,'color','black')
        xlim(opt.sta_corr_range_sec)
        xlabel('Corr time lag (red-green, sec)')
        ylabel('Counts')
        title(['Median. ' num2str(median(corr_lag_time),'%.2f') 'sec'])
        box off; axis square
        
        subplot(num_plot_rows,num_plot_cols,3+num_plot_cols*2); hold on
        histogram(peak_lag_time,'binwidth',.5,'facecolor','black','edgecolor','none')
        plot(median(peak_lag_time).*[1,1],ylim,'color','black')
        xlim(opt.sta_corr_range_sec)
        xlabel('Peak time lag (red-green, sec)')
        ylabel('Counts')
        title(['Median. ' num2str(median(peak_lag_time),'%.2f') 'sec'])
        box off; axis square
        
        this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_GreenRedCorr_' fd];
        suptitle(strrep(this_title,'_',' '))
        export_fig([this_save_path filesep this_title '.png'])
        
    end
    %% save correlation struct
    save([this_save_path filesep 'corr_struct'],'corr_struct')
end

if IF_GET_SPONTANEOUS
    %% --- spontaneous activity analysis ---
    % parames for event detection
    session_length = getOr(opt,'session_length',50);% minutes, for BE it's 50 minutes
    opt.event_thresh = 0;
    opt.baseline_frames_ed = round(0.5*opt.frame_rate); % sta baseline for event detection (ed)
    opt.pre_frames_ed = round(opt.frame_rate); % sta pre frames for event detection (ed)
    opt.post_frames_ed = 2*round(opt.frame_rate); % sta post frames  for event detection (ed)
    opt.IfPlot = 0; % plot event detection process
    opt.trig_peak_min = 0; % 90% - 10% rising signal level
    opt.fil_window = round(0.01*opt.frame_rate); % discard events that are too close to the preceding one
    opt.sensor_tau = .5; % sec, for event detection template
    opt.Th = 4; % for template matching threshold (4 was chosen in the Clements and Bekkers paper)
    opt.peri_event_window = 5; % sec. time window centered at behavioural events that will be cut out from trace for spontaneous activity analysis
    opt.min_spont_period = 10; %sec
    %% get 'spontaneous' activity periods
    % from MedPC start to session end
    % concatenate episodes with recorded behavioural events
    plot_color = color.(event_detection_channel);
    try
        input_trace = proc_traces.(['db_' event_detection_channel])(round(opt.medpc_start_time*opt.frame_rate+[1:session_length*opt.frame_rate*60])); % debleached GRAB-DA signal (zscored)
    catch
        disp('not enough frames')
    end
    peri_event_frames = opt.peri_event_window*opt.frame_rate;
    if ~isempty(event_frames)
        temp = struct2cell(event_frames);
        all_event_frames = unique(cell2mat(temp(~cellfun(@(x)isempty(x),temp)))); % frames with recorded behavioural events
        min_spont_period = round(opt.min_spont_period*opt.frame_rate); % minimum duration of a single spontaenous period
        [spont_trace] = get_spont_trace(input_trace,all_event_frames,peri_event_frames,min_spont_period);
    else
        spont_trace = input_trace;
    end
    %% event detection with template matching and iterative filtering
    % create template
    t = 0:1:opt.sensor_tau*opt.frame_rate*2;
    template = [0 exp(-t./(opt.sensor_tau.*opt.frame_rate))];
    smooth_window =  round(0.5*opt.frame_rate);
    smooth_func = @(x)movmean(x,smooth_window); % smoothing function
    num_iteration = 3; % number of iteration for event filtering
    smooth_trace = smooth_func(spont_trace);
    % downsample trace and template to speed up
    ds_factor = round(opt.frame_rate/50);%downsample to 50 Hz
    [ds_opt] = get_ds_opt(opt,ds_factor); % downsampled opt
    smooth_trace = downsample(smooth_trace,ds_factor);
    template = downsample(template,ds_factor);
    raw_events = template_matching(smooth_trace,template(:),ds_opt.Th);
    raw_event_trace = zeros(size(smooth_trace)); raw_event_trace(raw_events) = 1;
    [filt_events,eta_traces,event_amps] = filt_and_align_events(smooth_trace,raw_events,num_iteration,ds_opt);
    spont_activity_duration = numel(find(~isnan(spont_trace)))/opt.frame_rate/60; % minutes
    tot_recording_duration = length(input_trace)/opt.frame_rate/60; % minutes
    spont_event_rate = numel(filt_events)/(spont_activity_duration*60);% Hz
    %% plot detected events
    figure('name','spontaneous events','units','normalized','outerposition',[0 0 1 1])
    % show trace
    num_plot_rows = 2;num_plot_cols = 3;plot_count = 1;
    subplot(num_plot_rows,num_plot_cols,1:num_plot_cols); hold on;
    xticks = (1:length(smooth_trace))/ds_opt.frame_rate/60;
    plot(xticks,smooth_trace,'color',plot_color,'linewidth',1.5);
    % marked detected events
    arrayfun(@(x)plot([1 1].*(x+xticks(1)),ylim,'color',[.5 .5 .5]),filt_events./ds_opt.frame_rate/60)
    y_label = '\Delta F (zscored)'; xlabel('Time from session start (minutes)')
    title(['Spont. duration: ' num2str(spont_activity_duration,'%.2f') ' minutes; Spont. event rate: ' num2str(spont_event_rate,'%.2f'),' Hz'])
    % Zoomed in trace
    subplot(num_plot_rows,num_plot_cols,num_plot_cols+plot_count); hold on;
    plot_range = find(~isnan(smooth_trace)); plot_range = plot_range(1)+[1:60*ds_opt.frame_rate]; % show 1 minute after the onset of the first spontanoeus period
    xticks = (plot_range)/ds_opt.frame_rate;
    plot(xticks,smooth_trace(plot_range),'color',plot_color,'linewidth',1.5);
    % marked detected events
    arrayfun(@(x)plot([1 1].*(x+xticks(1)),ylim,'color',[.5 .5 .5]),(filt_events(filt_events>plot_range(1)&filt_events<plot_range(end))-plot_range(1))./ds_opt.frame_rate)
    y_label = '\Delta F (zscored)'; xlabel('Time from first spont. episode start(sec)')
    plot_count = plot_count+1;
    % STA
    subplot(num_plot_rows,num_plot_cols,num_plot_cols+plot_count); hold on;
    this_xticks = (-ds_opt.pre_frames_ed:1:0:1:ds_opt.post_frames_ed)./ds_opt.frame_rate;
    shadedErrorBar(this_xticks,nanmean(eta_traces,1), nanstd(eta_traces,[],1),{'color',plot_color,'linewidth',2},0.1);
    plot([0 0],ylim,'color',plot_color,'linewidth',2,'linestyle',':');
    box off; axis square;
    xlim([this_xticks(1) this_xticks(end)])
    xlabel('Time from event onset (sec)'); ylabel(y_label)
    title({[num2str(size(eta_traces,1)) ' events' ]})
    plot_count = plot_count+1;
    % Event amplitude
    subplot(num_plot_rows,num_plot_cols,num_plot_cols+plot_count); hold on;
    histogram(event_amps,'facecolor',[.5 .5 .5]); axis square
    xlabel('Event amplitude (\Delta F (zscored))')
    ylabel('Count')
    
    this_title = [this_mouse '_' this_date '_' session_type '_' this_location '_SpontaneousEvents_' event_detection_channel];
    suptitle(strrep(this_title,'_',' '))
    export_fig([this_save_path filesep this_title '.png'])
    
    %% save results
    spont_event_struct = struct();
    spont_event_struct.spont_trace = spont_trace;
    spont_event_struct.spont_activity_duration = spont_activity_duration;
    spont_event_struct.spont_event_rate = spont_event_rate;
    spont_event_struct.filt_events = filt_events;
    spont_event_struct.eta_traces = eta_traces;
    spont_event_struct.event_amps = event_amps;
    spont_event_struct.ds_opt = ds_opt;
    save([this_save_path filesep 'spont_event_struct'],'spont_event_struct')
    
end


%% --- close figures ---
close all
end

