function [ event_frames,this_tta_raw,all_event_amp] = filt_and_align_events( interp_trace,event_frames,num_iteration,options)
% align and filter events (modified from Plasticity/filt_and_align_tta)
% interp_traces: Fcell - Fneu, with photostim frames interpolated
% event_frames : frame indices of the thresholded suite2p spike train (thresh = mean + 2 std)
% num_iteration: 3 works fine
pre_frames               = options.pre_frames_ed;
post_frames              = options.post_frames_ed;
baseline_frames          = [-options.baseline_frames_ed:1:0]+options.pre_frames_ed+1; % normalise wrt these frames to make stas
search_peak_frames       = options.pre_frames_ed; % search peak in sta
IfNormalise              = getOr(options,'IfNormalise',false);
% trig_peak_min            = options.trig_peak_min;
trig_peak_min            = getOr(options,'trig_peak_min',0); % changed to signal level difference 90% - 10%  ZZ 2022 
IfPlot                   = getOr(options,'IfPlot',0);
fil_window               = getOr(options,'fil_window',5); 
event_frames = sort(event_frames);
% define average filter
windowSize = 5;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
org_event = event_frames;
org_trace = interp_trace;

discarded_traces = [];
all_event_amp = [];
dis_idx = 1;

% check for debugging
if(IfPlot)
    figure('name','event detection','units','normalized','outerposition',[0 0 .5 1])
    if IfNormalise
        base_dff = nanmean(interp_trace(interp_trace<quantile(interp_trace,0.30)));
        dff = (interp_trace-base_dff)./base_dff;
    else
        dff = interp_trace;
    end
end
for i = 1:num_iteration
    if isempty(event_frames)
        all_event_amp = [];
        break
    end
    event_frames = unique(event_frames);
%     [~,~,~,~,~,raw_sta]...
%         = make_sta_from_traces( interp_trace,event_frames,pre_frames,post_frames,baseline_frames);
%     this_tta_raw = squeeze(raw_sta);
%      
    % check event detection results in each iteration
    if(IfPlot)
        subplot(num_iteration+1,1,i)
        hold on;
        plot(dff,'color','black');
        spike_train = zeros(size(dff));
        spike_train(event_frames) = max(dff);
        plot(spike_train,'color','r');
        xlim([1 length(interp_trace)]);
        ylabel(['Iteration' num2str(i)]);
    end
    
    % discard events too close to the preceding one
    temp_frames = event_frames;
    for event_idx = 2:numel(event_frames)
        if(temp_frames(event_idx)-temp_frames(event_idx-1)<fil_window)
            temp_frames(event_idx) = nan;
%             discarded_traces = [discarded_traces;this_tta_raw(event_idx,:)];
        end
    end
    event_frames = temp_frames(~isnan(temp_frames));
    
    [~,~,~,~,~,raw_sta]...
        = make_sta_from_traces( interp_trace,event_frames,pre_frames,post_frames,baseline_frames);
    this_tta_raw = squeeze(raw_sta);

    % Filter out small events
    if(~isempty(this_tta_raw))
        if(size(this_tta_raw,2)~=1)
            this_indices = [];
            this_onset = [];
            for m = 1:size(this_tta_raw,1)
                this_trace = this_tta_raw(m,1:pre_frames+search_peak_frames);
                if(any(isnan(this_trace)))
                    continue
                end
                this_trace = filtfilt(b,a,medfilt1(this_trace));     % smooth event traces before getting rise time
                [this_rise_time,lwrCross,uprCross,lwrRef,uprRef] = risetime(this_trace,'PercentReferenceLevels',[10 90]);
                if(isempty(this_rise_time))
                    discarded_traces = [discarded_traces;this_tta_raw(m,:)];
                    continue
                else
                    [~,temp_idx] = max(this_rise_time);
                end
                if((uprRef-lwrRef)<trig_peak_min) 
                    discarded_traces = [discarded_traces;this_tta_raw(m,:)];
                    dis_idx = dis_idx+1;
                    continue
                else
                    this_indices = [this_indices,m];
                    this_onset = [this_onset,round(lwrCross(temp_idx))];
                    temp = event_frames(m);
                    event_frames(m) = temp+round(lwrCross(temp_idx))-pre_frames-1;
                end
                
            end
            event_frames = event_frames(this_indices);
            
            
            
            % discard events too close to the preceding one
            temp_frames = event_frames;
            for event_idx = 2:numel(event_frames)
                if(temp_frames(event_idx)-temp_frames(event_idx-1)<fil_window)
                    temp_frames(event_idx) = nan;
                    %             discarded_traces = [discarded_traces;this_tta_raw(event_idx,:)];
                end
            end
            event_frames = temp_frames(~isnan(temp_frames));
            
        else
            this_tta_raw = this_tta_raw';
        end
        if i == 1
            org_traces = this_tta_raw;
        end
    else
        event_frames = [];
        return
    end
end

% adjust onset time for the last time
if(~isempty(event_frames))
    [~,~,~,~,~,raw_sta]...
        = make_sta_from_traces( interp_trace,event_frames,pre_frames,post_frames,baseline_frames);
    this_tta_raw = squeeze(raw_sta);
end
% adjust event frame to the 20% rise time
if(~isempty(this_tta_raw))
    all_event_amp = zeros(1,size(this_tta_raw,2));
    if(size(this_tta_raw,2)~=1)
        this_indices = [];
        this_onset = [];
        for m = 1:size(this_tta_raw,1)
            this_trace = this_tta_raw(m,1:pre_frames+search_peak_frames);
            if(any(isnan(this_trace)))
                continue
            end
            this_trace = filtfilt(b,a,medfilt1(this_trace));     % smooth event traces before getting rise time
            [this_rise_time,lwrCross,uprCross,lwrRef,uprRef] = risetime(this_trace,'PercentReferenceLevels',[10 90]);
            if(isempty(this_rise_time))
                discarded_traces = [discarded_traces;this_tta_raw(m,:)];
                continue
            else
                [~,temp_idx] = max(this_rise_time);
            end
            if((uprRef-lwrRef)<trig_peak_min)
                discarded_traces = [discarded_traces;this_tta_raw(m,:)];
                dis_idx = dis_idx+1;
                continue
            else
                this_indices = [this_indices,m];
                this_onset = [this_onset,round(lwrCross(temp_idx))];
                temp = event_frames(m);
                event_frames(m) = temp+round(lwrCross(temp_idx))-pre_frames-1; % was lwrCross,uprCross
%                 event_frames(m) = temp+round(uprCross(temp_idx))-pre_frames-3; % was lwrCross,uprCross
%                 event_frames(m) = temp+max_frame-pre_frames-3; % was lwrCross,uprCross
                all_event_amp(m)= uprRef;
            end            
        end
        event_frames = event_frames(this_indices);
        all_event_amp = all_event_amp(this_indices);
    else
        this_tta_raw = this_tta_raw';
    end    
    event_frames = unique(event_frames);
    event_frames(event_frames<0) = [];
    % discard events too close to the preceding one
    temp_frames = event_frames;
    for event_idx = 2:numel(event_frames)
        if(temp_frames(event_idx)-temp_frames(event_idx-1)<fil_window)
            temp_frames(event_idx) = nan;
%             discarded_traces = [discarded_traces;this_tta_raw(event_idx,:)];
        end
    end
    event_frames = temp_frames(~isnan(temp_frames));

    [~,~,~,~,~,raw_sta]...
        = make_sta_from_traces( interp_trace,event_frames,pre_frames,post_frames,baseline_frames);
    this_tta_raw = squeeze(raw_sta);
    
    if(IfPlot)
        subplot(num_iteration+1,1,num_iteration+1)
        hold on;
        plot(dff,'color','black');
        spike_train = zeros(size(dff));
        spike_train(event_frames) = max(dff);
        plot(spike_train,'color','r');
        xlim([1 length(interp_trace)]);
        ylabel('Final result');
        suptitle(['Event threshold = ' num2str(trig_peak_min)])
    end
    
else
    event_frames = [];
    all_event_amp = [];
    return
end


if(IfPlot)
    ylimit(1) = min(min([org_traces;org_traces]));
    ylimit(2) = max(max([org_traces;org_traces]));
    xaxis_value = (-pre_frames:post_frames);
    % plot check
    figure('name','aligning trigger events')
    subplot(2,3,1)
    title('before filtering')
    hold on
    imagesc(org_traces,ylimit)
    plot([pre_frames pre_frames],ylim,'color','r');
    xlim(size(org_traces(1,:)));
    
    subplot(2,3,2)
    title('after filtering')
    if(~isempty(this_tta_raw))
        hold on
        imagesc(this_tta_raw,ylimit)
        plot([pre_frames pre_frames],ylim,'color','r');
        xlim(size(org_traces(1,:)));
        colormap gray
    end
    
    subplot(2,3,3)
    title('discarded small events')
    if(~isempty(discarded_traces))
        hold on
        for i = 1:size(discarded_traces,1)
            plot(xaxis_value,discarded_traces(i,:)+i,'color',[0.3 0.3 0.3]);
        end
        xlim([xaxis_value(1) xaxis_value(end)])
        plot([0 0],ylim,'color','r');
    end
    suptitle(['dF/F threshold = ' num2str(trig_peak_min)])
    
%     figure('name','average trace')
    subplot(2,3,4)
    title(['before filtering: ',num2str(size(org_traces,1))])
    hold on
    for i = 1:size(org_traces,1)
        plot(xaxis_value,org_traces(i,:));
    end
    plot(xaxis_value,mean(org_traces,1),'color','black','linewidth',2)
    %     shadedErrorBar(xaxis_value,mean(org_traces,1),nanstd(org_traces,0,1)./sqrt(size(org_traces,1)),{'color',[0.3 0.3 0.3],'markerfacecolor',[0.3 0.3 0.3]},0.5);
    xlim([xaxis_value(1) xaxis_value(end)])
    ylim(ylimit);
    axis square
    
    subplot(2,3,5)
    title(['after filtering: ', num2str(size(this_tta_raw,1))])
    hold on
    for i = 1:size(this_tta_raw,1)
        plot(xaxis_value,this_tta_raw(i,:));
    end
    plot(xaxis_value,mean(this_tta_raw,1),'color','black','linewidth',2)
    %     shadedErrorBar(xaxis_value,mean(this_tta_raw,1),nanstd(this_tta_raw,0,1)./sqrt(size(this_tta_raw,1)),{'color',[0.3 0.3 0.3],'markerfacecolor',[0.3 0.3 0.3]},0.5);
    xlim([xaxis_value(1) xaxis_value(end)])
    ylim(ylimit);
    axis square
    

    suptitle(['Event threshold = ' num2str(trig_peak_min)])
    
    
end
end

