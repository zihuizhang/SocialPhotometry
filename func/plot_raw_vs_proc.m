function [] = plot_raw_vs_proc(this_red,this_db_red,plot_color,proc_type,opt, varargin)
y_label = 'Zscored';
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'y_label')
        y_label = varargin{v+1};
    end
end
ds_factor = round(opt.frame_rate/50); % downsample to 50Hz
figure('name','raw vs processed data','units','normalized','outerposition',[0 0 1 .5]); 
subplot(2,9,1:9);hold on;
h(2) = plot(downsample(medfilt1(this_db_red,10), ds_factor),'color',plot_color,'displayname',proc_type);
h(1) = plot(downsample(medfilt1(this_red,10), ds_factor),'color',[.5 .5 .5],'displayname','Raw');
xlim([round(opt.discard_first_frames/ds_factor), round(opt.num_frames/ds_factor)]);
ylabel(y_label); legend
title(['Full session, ' num2str(round(opt.num_frames/opt.frame_rate)), ' sec, sample rate ' num2str(round(opt.frame_rate)) ' Hz, downsampled by ' num2str(ds_factor) ,'x'])

subplot(2,9,10:12);hold on; % zoomed in first 15 seconds
h(1) = plot(this_red,'color',[.5 .5 .5],'displayname',proc_type);
h(2) = plot(this_db_red,'color',plot_color,'displayname','Debleached');
xlim([opt.discard_first_frames, round(15*opt.frame_rate)]);
ylabel(y_label)
title('First 15 sec')

subplot(2,9,13:15);hold on; % zoomed in first 15 seconds
h(1) = plot(this_red,'color',[.5 .5 .5],'displayname',proc_type);
h(2) = plot(this_db_red,'color',plot_color,'displayname','Debleached');
xlim([round(opt.num_frames/2-opt.frame_rate*7.5) round(opt.num_frames/2+opt.frame_rate*7.5)]);
title('Mid 15 sec')

subplot(2,9,16:18);hold on; % zoomed in last 15 seconds
h(1) = plot(this_red,'color',[.5 .5 .5],'displayname',proc_type);
h(2) = plot(this_db_red,'color',plot_color,'displayname','Debleached');
xlim([round(opt.num_frames-15*opt.frame_rate) opt.num_frames]);
title('Last 15 sec')
end

