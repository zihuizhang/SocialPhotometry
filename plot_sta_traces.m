function [] = plot_sta_traces(sta_struct,plot_trace_fd,plot_channels,plot_fds,ylimits,y_label,color,opt,varargin)
IF_FILT = false;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_FILT')
        IF_FILT = varargin{v+1};
    end
end
filt_size = 3;
num_plot_cols = numel(plot_fds); num_plot_rows = numel(plot_channels)+1; 
this_xticks = (-opt.pre_frames:1:0:1:opt.post_frames)./opt.frame_rate;
photo_duration_sec = getOr(opt,'photo_duration_sec',[]);
for f = 1:num_plot_cols
    % shaded error bar traces
    subplot(num_plot_rows,num_plot_cols,f); hold on
    for c = 1:numel(plot_channels)
        this_channel = plot_channels{c}; fd = plot_fds{f};
        this_traces = sta_struct.([this_channel plot_trace_fd]).(fd);
        if IF_FILT
            for t = 1:size(this_traces,1)
                this_traces(t,:) = movmean( medfilt1(this_traces(t,:),filt_size),filt_size);
            end
        end
        shadedErrorBar(this_xticks,nanmean(this_traces,1), nanstd(this_traces,[],1)./sqrt(size(this_traces,1)),{'color',color.(this_channel),'linewidth',2},0.1);
   
    end
%     ylim([-2 5]) % for example plots
    plot([0 0],ylim,'color',color.(fd),'linewidth',2,'linestyle',':'); 
    if ~isempty(photo_duration_sec) % photostim off time
        plot([1 1].*photo_duration_sec,ylim,'color',[.5 .5 .5],'linewidth',1.5,'linestyle',':');
    end
    box off; axis square;
    xlim([-opt.sta_pre_sec opt.sta_post_sec])
        xlabel('Time from event (sec)'); ylabel(y_label)
    title({strrep(fd,'_',' '); [num2str(size(this_traces,1)) ' trials']})   

    
    % imagesc
    cc = 1;
    for c = 1:numel(plot_channels)
        this_channel = plot_channels{c};
%         if strcmp(this_channel,'uv') % don't bother uv channel
%             continue
%         end
        
        subplot(num_plot_rows,num_plot_cols,f+num_plot_cols*cc); hold on
        this_traces = sta_struct.([this_channel plot_trace_fd]).(fd);
        if ~isempty(this_traces)
        if isstruct(ylimits)
            ylimit = ylimits.(this_channel);
        else
            ylimit = ylimits;
        end
        imagesc(this_traces);colormap(b2r(ylimit(1),ylimit(2)));% axis square
        plot([1 1].*(opt.pre_frames+1),xlim,'color',[.5 .5 .5])
        if ~isempty(photo_duration_sec) % photostim off time
            plot([1 1].*(photo_duration_sec.*opt.frame_rate+opt.pre_frames),ylim,'color',[.5 .5 .5],'linewidth',1.5,'linestyle',':');
        end
        colorbar('eastoutside')
        ylim([0.5 size(this_traces,1)+0.5]); xlim([1, size(this_traces,2)])
        ylabel('Trials'); title(strrep(this_channel,'_',' '))
        set(gca,'xtick',[],'xcolor','w')
        else
            axis off
        end
        cc = cc+1;
    end
       
end

end

