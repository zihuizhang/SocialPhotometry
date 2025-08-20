function [ P,x_positions,y_pos,stats] = scatter_cmp_conditions(values,plot_name,subplot,input_color_lut,varargin)
% scatter plot of fields in values
% inputs: values,plot_name,subplot,input_color_lut,varargin
% value.fields are row vectors
% TO DO: compare more than 2 fields
flag_color = 0; markerfacealpha = 1;
xtick_label = [];
xlimit = [];ylimit = []; connect_pairs = [];
if(nargin>2)
    ifsubplot = subplot;
    if(~ifsubplot)
        figure('Name',plot_name)
    end
    try
        temp = input_color_lut;
        if(~isempty(temp))
            flag_color = 1;
            color_lut = input_color_lut;
        else
            color_lut = repmat([0 0 0],numel(fields(values)),1);
        end
    catch
        color_lut = [0.3 0.3 0.3;[231 73 24]./255;0.3 0.3 0.3];
    end
else
    figure('Name',plot_name)
end

par = inputParser;
defaultIfPrePostOnly = 0;
default_test_type = 'ranksum';
addParameter(par,'connect_scatter',1);
addParameter(par,'PlotPrePostOnly',defaultIfPrePostOnly);
addParameter(par,'test_type',default_test_type);
addParameter(par,'IfAvgByAnimal',0);
addParameter(par,'IfBoxWhisker',0);
addParameter(par,'IfViolin',0);
addParameter(par,'x_shift',0);
addParameter(par,'x_interval',1);
addParameter(par,'barwidth',.5);
addParameter(par,'barEdgecolor',[0 0 0]);
addParameter(par,'add_jitter',0);
addParameter(par,'add_yjitter',0);
addParameter(par,'animal_colors',[]);
addParameter(par,'IfColorAnimals',0);
addParameter(par,'tail','both');
addParameter(par,'IfHistogram',0);
addParameter(par,'plot_stats',0);
addParameter(par,'BinWidth',10);
addParameter(par,'displayopt','off');
addParameter(par,'Normalization','probability')
addParameter(par,'DisplayStyle','stairs')
addParameter(par,'BriefXlabel',0)
addParameter(par,'ShowMeanInXlabel',0)
addParameter(par,'VeryBriefXlabel',0)
addParameter(par,'xtick_label',[])
addParameter(par,'xlimit',[])
addParameter(par,'ylimit',[])
addParameter(par,'plot_mean',1) % for violin plot
addParameter(par,'plot_median',1) % for violin plot
addParameter(par,'barfacecolor','none')
addParameter(par,'markerfacealpha',1)
addParameter(par,'plot_se',1)
addParameter(par,'plot_scatter',1)
addParameter(par,'connect_pairs',[])
addParameter(par,'factors',[]) % for 2-way anova; field names format: factor1_factor2

parse(par,varargin{:})
if(~ isempty(par.Results.factors))
    factors = par.Results.factors;
end
if(~ isempty(par.Results.xlimit))
    xlimit = par.Results.xlimit;
end
if(~ isempty(par.Results.ylimit))
    ylimit = par.Results.ylimit;
end
if(~ isempty(par.Results.tail))
    tail = par.Results.tail;
end
if(~isempty(par.Results.PlotPrePostOnly))
    IfPlotPrePostOnly = par.Results.PlotPrePostOnly;
end
if(~isempty(par.Results.test_type))
    test_type = par.Results.test_type;
end
if(~isempty(par.Results.connect_scatter))
    connect_scatter = par.Results.connect_scatter;
end
if(~isempty(par.Results.connect_pairs))
    connect_pairs = par.Results.connect_pairs;
end
if(~isempty(par.Results.IfBoxWhisker))
    IfBoxWhisker = par.Results.IfBoxWhisker;
end
if(~isempty(par.Results.x_shift))
    x_shift = par.Results.x_shift;
end
if(~isempty(par.Results.x_interval))
    x_interval = par.Results.x_interval;
end
if(~isempty(par.Results.barwidth))
    barwidth = par.Results.barwidth;
end
if(~isempty(par.Results.barEdgecolor))
    barEdgecolor = par.Results.barEdgecolor;
end
if(~isempty(par.Results.add_jitter))
    IfAddJitter = par.Results.add_jitter;
end
if(~isempty(par.Results.add_jitter))
    IfAddYJitter = par.Results.add_yjitter;
end
if(~isempty(par.Results.animal_colors))
    animal_colors = par.Results.animal_colors;
end

if(~isempty(par.Results.IfColorAnimals))
    IfColorAnimals = par.Results.IfColorAnimals;
end

if(~isempty(par.Results.IfHistogram))
    IfHistogram = par.Results.IfHistogram;
end

if(~isempty(par.Results.plot_stats))
    plot_stats = par.Results.plot_stats;
end

if(~isempty(par.Results.BinWidth))
    BinWidth = par.Results.BinWidth;
end

if(~isempty(par.Results.displayopt))
    displayopt = par.Results.displayopt;
end

if(~isempty(par.Results.Normalization))
    Normalization = par.Results.Normalization;
end

if(~isempty(par.Results.DisplayStyle))
    DisplayStyle = par.Results.DisplayStyle;
end

if(~isempty(par.Results.BriefXlabel))
    BriefXlabel = par.Results.BriefXlabel;
end

if(~isempty(par.Results.BriefXlabel))
    ShowMeanInXlabel = par.Results.ShowMeanInXlabel;
end

if(~isempty(par.Results.VeryBriefXlabel))
    VeryBriefXlabel = par.Results.VeryBriefXlabel;
end

if(~isempty(par.Results.xtick_label))
    xtick_label = par.Results.xtick_label;
end

if(~isempty(par.Results.IfViolin))
    IfViolin = par.Results.IfViolin;
end
if(~isempty(par.Results.plot_mean))
    plot_mean = par.Results.plot_mean;
end
if(~isempty(par.Results.plot_median))
    plot_median = par.Results.plot_median;
end
if(~isempty(par.Results.barfacecolor))
    barfacecolor = par.Results.barfacecolor;
end
if(~isempty(par.Results.markerfacealpha))
    markerfacealpha = par.Results.markerfacealpha;
end
if(~isempty(par.Results.plot_se))
    plot_se = par.Results.plot_se;
end
if(~isempty(par.Results.plot_scatter))
    plot_scatter = par.Results.plot_scatter;
end
% plot each field in value
hold on
fd_names = fieldnames(values);
num_plot = numel(fd_names);

% initialise
max_size = 0;
for m = 1:num_plot
    temp_value = values.(fd_names{m});
    this_size = numel(temp_value);
    if(this_size>max_size)
        max_size = this_size;
    end
end
mat_all_values = nan(max_size,num_plot);
num_samples = zeros(1,num_plot);
num_nonnan_samples = zeros(1,num_plot);

for m =1:num_plot
    temp_value = values.(fd_names{m});
    if(size(temp_value,2)>1)
        temp_value = temp_value';   % rotate to row vector
    end
    %     animal_value = arrayfun(@(x)nanmean(extractfield(values(x),fields{m})),1:size(values,2));
    mat_all_values(1:length(temp_value),m)= temp_value;
    num_samples(m) = length(temp_value);
    num_nonnan_samples(m) = length(find(~isnan(temp_value)>0));
end

%% exclude rows that contain nans 
% rows_with_nans = arrayfun(@(x)sum(isnan(mat_all_values(x,:)))>0,1:size(mat_all_values,1));
% mat_all_values(rows_with_nans,:) = [];

% % stat stest on all combinations of fields
if num_plot>1
    all_cmp = nchoosek(1:num_plot,2);
else
    all_cmp =[1 1];
end
stats = struct();
max_value = max(mat_all_values(:));
min_value = min(mat_all_values(:));
if IfHistogram
    max_value = 0.5;
    min_value = 0;
end
if ~isempty(ylimit)
    max_value = ylimit(2);
    min_value = ylimit(1);
end
this_plot_stats_hight = (max_value-min_value)*0.1+max_value*0.9;
plot_stats_step = 0.2*max_value;


if(~(contains(test_type,'anova')||(strcmp(test_type,'KW')||(strcmp(test_type,'FM')))))
    % compare pairs - not corrected for multi comparison
    if ~isempty(connect_pairs)&&strcmp(test_type,'signrank')
        all_cmp = cell2mat(connect_pairs');
    end
    for c = 1:size(all_cmp,1)
        idx1 = all_cmp(c,1);
        idx2 = all_cmp(c,2);
        stats(c).field1 = fd_names{idx1};
        stats(c).field2 = fd_names(idx2);
        if(isempty(test_type))
            [P,h] = ranksum(mat_all_values(:,idx1),mat_all_values(:,idx2));
        else
            switch test_type
                case 'ttest'
                    [h,P]=ttest(mat_all_values(:,idx1),mat_all_values(:,idx2));
                case 'ranksum'
                    try
                        [P,h] = ranksum(mat_all_values(:,idx1),mat_all_values(:,idx2));
                    catch
                        P = 1; h = 0;
                    end
                case 'ttest2'
                    [h,P]=ttest2(mat_all_values(:,idx1),mat_all_values(:,idx2),'tail',tail,'Vartype','unequal');
                case 'signrank'
                    try
                        [P,h] = signrank(mat_all_values(:,idx1),mat_all_values(:,idx2),'tail',tail);
                    catch
                        P = 1; h = 0;
                    end
                case 'KS'
                    [h,P] = kstest2(mat_all_values(:,idx1),mat_all_values(:,idx2));
                    %                 IfHistogram = 1; % will plot histogram if chosen KS test
            end
        end
        stats(c).P = P;
        stats(c).h = h;

        % plot asterisks
        if(plot_stats&&(~isempty(plot_stats_step)))
            if P>0.05 % only show significant
                continue
            end
            this_plot_stats_hight = this_plot_stats_hight+ plot_stats_step;
            plot([idx1 idx2],this_plot_stats_hight.*[1,1],'black');
            text(mean([idx1 idx2]),this_plot_stats_hight, get_stars_for_P(P), 'VerticalAlignment','bottom','HorizontalAlignment','center');
        end

    end
else
    P = [];
    % compare pairs - not corrected for multi comparison
    [P,tbl,mult_stats] = multcmp_stats_on_struct( values, 'test_type',test_type,'displayopt',displayopt);
    if ~(strcmp(test_type,'anovan')||strcmp(test_type,'FM'))
        [C,means,h,gnames] = multcompare(mult_stats,'Display',displayopt);
    else
        [C,means,h,gnames] = multcompare(mult_stats,'Display',displayopt,'Dimension',[1 2]);
    end
    for c = 1:size(C,1)
        stats(c).P = C(c,end);
        idx1 = C(c,1);
        idx2 = C(c,2);
        fd_name1 = gnames{idx1}; fd_name1 = strrep(fd_name1,'X1=','');fd_name1 = strrep(fd_name1,',X2=','_');
        fd_name2 = gnames{idx2}; fd_name2 = strrep(fd_name2,'X1=','');fd_name2 = strrep(fd_name2,',X2=','_');
        idx1 = find(cellfun(@(x)strcmp(x,fd_name1),fd_names));
        idx2 = find(cellfun(@(x)strcmp(x,fd_name2),fd_names));
        stats(c).pair_idx = [idx1,idx2];

        % plot asterisks
        if(plot_stats)&&(stats(c).P<0.05)
            this_plot_stats_hight = this_plot_stats_hight+ plot_stats_step;
            plot([idx1 idx2],this_plot_stats_hight.*[1,1],'black');
            text(mean([idx1 idx2]),this_plot_stats_hight, get_stars_for_P( stats(c).P), 'VerticalAlignment','bottom','HorizontalAlignment','center');

        end
    end
end



try
    if(isempty(test_type))
        [P,h] = ranksum(mat_all_values(:,1),mat_all_values(:,end));
    else
        switch test_type
            case 'ttest'
                [h,P]=ttest(mat_all_values(:,1),mat_all_values(:,end));
            case 'ranksum'
                [P,h] = ranksum(mat_all_values(:,1),mat_all_values(:,end));
            case 'ttest2'
                [h,P]=ttest2(mat_all_values(:,1),mat_all_values(:,end));
            case 'signrank'
                [P,h] = signrank(mat_all_values(:,1),mat_all_values(:,end),'tail',tail);
        end
    end
catch
    P = 1; h = 0;
    warning('stats test error')
end

mean_values = nanmean(mat_all_values,1);
diff_value = 100.*(mean_values(:,end)-mean_values(:,1))./abs(mean_values(:,1));
se_values = nanstd( mat_all_values,0,1)./sqrt(arrayfun(@(x)numel(find(~isnan(mat_all_values(:,x)))),1:size(mat_all_values,2)));
sd_values = nanstd( mat_all_values,0,1);
sk_values = skewness(mat_all_values);
median_values = nanmedian(mat_all_values,1);
%% plot
x_positions = (1:num_plot).*x_interval+x_shift;
if ~IfPlotPrePostOnly   % plot all conditions
    if((~IfBoxWhisker)&&(~IfHistogram)) % mean +- se
        hold on
        if(connect_scatter)&&(~IfColorAnimals) % connect dots with gray lines
            if isempty(connect_pairs)
                plot(repmat(x_positions,size(mat_all_values,1),1)',mat_all_values','color',[.7 .7 .7],'linewidth',0.5)
            else
                for c = 1:numel(connect_pairs)
                    plot(repmat(x_positions(connect_pairs{c}),size(mat_all_values,1),1)',mat_all_values(:,connect_pairs{c})','color',[.7 .7 .7],'linewidth',0.5)
                end

            end
        end
        if IfColorAnimals % connect dots colored by animal
            temp_value = mat_all_values(:,i);
            for j = 1:numel(temp_value)
                if isempty(connect_pairs)
                    plot(x_positions,mat_all_values(j,:),'color',animal_colors(j,:),'linewidth',0.25)
                else
                    for c = 1:numel(connect_pairs)
                        plot(x_positions(connect_pairs{c}),mat_all_values(j,connect_pairs{c}),'color',animal_colors(j,:),'linewidth',0.5)
                    end

                end
            end
        end
        for i =1:num_plot
            temp_value = mat_all_values(:,i);
            temp_value = temp_value(:);
            if isempty(temp_value)
                continue;
            end
            if IfViolin

                violin(temp_value,'x',x_positions(i).*ones(size(temp_value)),...
                    'plot_median',plot_median,'plot_mean',plot_mean,...
                    'facealpha',.3,'facecolor',color_lut(i,:),'edgecolor','none','mc',color_lut(i,:),'medc',color_lut(i,:),'LineWidth',1);

            else
                hold on
                % scatter individual dots
                if plot_scatter
                    if(~IfColorAnimals)
                        plotx = x_positions(i).*ones(size(temp_value));
                        ploty = temp_value;
                        if(IfAddJitter)
                            jitters = randi([-200 200],size(temp_value))./200.*x_interval.*0.1;
                            plotx = x_positions(i).*ones(size(temp_value))+jitters;
                        end
                        if IfAddYJitter
                            ploty = temp_value + randi([-200 200],size(temp_value))./200.*max(temp_value(:)).*0.05;
                        end
                        scatter(plotx,ploty,15,'MarkerFaceColor',color_lut(i,:),'MarkerEdgeColor','none','markerfacealpha',markerfacealpha);
                    else

                        for j = 1:numel(temp_value)
                            jitters = 0;
                            if IfAddJitter
                                jitters = randi([-200 200],1)./200.*x_interval.*0.1;
                            end
                            scatter(x_positions(i)+jitters,temp_value(j),15,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',animal_colors(j,:))
                        end
                    end
                end
                % bar and whisker
                if ~plot_se
                    barwitherr(sd_values(i),mean_values(i),'XData',x_positions(i),'barwidth',barwidth,'FaceColor',barfacecolor,'EdgeColor',color_lut(i,:),'LineWidth',1.5);
                else
                    barwitherr(se_values(i),mean_values(i),'XData',x_positions(i),'barwidth',barwidth,'FaceColor',barfacecolor,'EdgeColor',color_lut(i,:),'LineWidth',1.5);
                end
            end
        end
        if isempty(xtick_label)
            field_names = cellfun(@(x)strrep(x,'_',' '),fd_names,'UniformOutput',false);
        else
            field_names = xtick_label;
        end
        set(gca,'XTick',[x_positions(1)-x_interval, x_positions, x_positions(end)+x_interval],'XTickLabel',{'' char(field_names) ''});
    elseif IfBoxWhisker % plot interquantile range
        boxplot(mat_all_values,'Whisker',1,'colors',color_lut)
        if isempty(xtick_label)
            field_names = cellfun(@(x)strrep(x,'_',' '),fd_names,'UniformOutput',false);
        else
            field_names = xtick_label;
        end
        set(gca,'XTick',[x_positions(1)-x_interval, x_positions, x_positions(end)+x_interval],'XTickLabel',{'' char(field_names) ''});
    elseif IfHistogram
        hold on
        for i =1:num_plot
            temp_value = values.(fd_names{i});
            if(strcmp(DisplayStyle,'stairs'))
                hp(i) = histogram(temp_value,'Normalization',Normalization,'DisplayStyle',DisplayStyle,'facecolor','none', 'edgecolor',color_lut(i,:),'linewidth',1,'Binwidth',BinWidth);
            else
                hp(i) = histogram(temp_value,'Normalization',Normalization,'DisplayStyle',DisplayStyle,'facecolor',color_lut(i,:), 'edgecolor','none','Binwidth',BinWidth,'FaceAlpha',1);
            end
            plot([1,1].*nanmean(temp_value),ylim,'color',color_lut(i,:))
        end
        legend_name = {};
        for ii = 1:numel(fd_names)
            legend_name{ii} = strrep(fd_names{ii},'_', ' ');
        end
        if ~isempty(xtick_label)
            legend_name = xtick_label;
        end
        if~isempty(xlimit)
            xlim(xlimit)
        end
        legend(hp(1:num_plot),legend_name,'Location','northoutside','box','off')
    end



else  % plot only the first and last set of values
    barwitherr(se_values(1),mean_values(1),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',.5)
    barwitherr(se_values(end),mean_values(end),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',.5)

    num_plot = 2;
    plot([mat_all_values(:,1), mat_all_values(:,3)]','color',[.7 .7 .7],'linewidth',0.5)
    temp_value = values.pre';
    scatter(ones(1,length(temp_value)),temp_value,50,'MarkerFaceColor',color_lut(1,:),'MarkerEdgeColor',color_lut(i,:),...
        'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    temp_value = values.post';
    scatter(2.*ones(1,length(temp_value)),temp_value,50,'MarkerFaceColor',color_lut(2,:),'MarkerEdgeColor',color_lut(i,:),...
        'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    set(gca,'XTick',0:3,'XTickLabel',{'' 'pre' 'post' ''});

end

if(~IfHistogram)
    ylabel(plot_name)
    xlim([0 num_plot+1])
else
    plot_name = strrep(plot_name,'_',' ');
    title(plot_name)
    ylabel(Normalization)
end
% xtickangle(45)
% show stats if only two fields are compared
if ~ BriefXlabel
    if(num_plot == 2)
        xlabel({[ test_type ' P=' num2str(P,'%10.1e') ' h= ' num2str(h)];...
            ['#Samples: ' num2str(num_samples)];...
            ['#Numeric: ' num2str(num_nonnan_samples)];...
            ['Mean: '  num2str(mean_values(1),'%10.3f') '; '  num2str(mean_values(end),'%10.3f')];...
            ['Median:' num2str(median_values(1),'%10.3f') '; '  num2str(median_values(end),'%10.3f')]
            ['SE: '  num2str(se_values(1),'%10.3f') '; '  num2str(se_values(end),'%10.3f') ];...
            ['SD: '  num2str(sd_values(1),'%10.3f') '; '  num2str(sd_values(end),'%10.3f') ];...
            ['Skewness: '  num2str(sk_values(1),'%10.3f') '; '  num2str(sk_values(end),'%10.3f')];...
            ['%Difference: ' num2str(diff_value,'%10.3f')]})
    else
        xlabel({
            ['#Samples: ' num2str(num_samples)];...
            ['#Numeric: ' num2str(num_nonnan_samples)];...
            ['Mean:  ' num2str(mean_values,'%10.3f')];...
            ['Median:' num2str(median_values,'%10.3f')];...
            ['SE:    '  num2str(se_values,'%10.3f')];...
            ['SD: '  num2str(sd_values,'%10.3f')];...
            ['Skew:  '  num2str(sk_values,'%10.3f')] })

    end
elseif ShowMeanInXlabel

    xlabel({['#Samples: ' num2str(num_samples)];...
        ['#Numeric: ' num2str(num_nonnan_samples)];...
        ['Mean: '  num2str(mean_values,'%10.2f') ];...
        ['SD: '  num2str(sd_values,'%10.2f')]})
else
    if (num_plot == 2)
        xlabel({ [ test_type ' P=' num2str(P,'%10.1e') ' h= ' num2str(h)];...
            ['#Samples: ' num2str(num_nonnan_samples)]
            ['Mean: '  num2str(mean_values,'%10.3f') ];...
            ['SD: '  num2str(sd_values,'%10.3f')]})

    else

        xlabel({
            ['#Numeric: ' num2str(num_nonnan_samples)];...
            ['Mean: '  num2str(mean_values,'%10.3f') ];...
            ['SE:    '  num2str(se_values,'%10.3f')];...
            ['SD: '  num2str(sd_values,'%10.3f')]})
    end

    %     xlabel({[ test_type ' P=' num2str(P,'%10.1e') ' h= ' num2str(h)]})
end
if VeryBriefXlabel
    xlabel({['Mean: '  num2str(mean_values,'%10.3f')]})
end
set(gcf,'color','w')
set(gca,'TickDir','out');
y_pos = double(1.1*max(max(mat_all_values)));
end

