function [this_filt_traces] = filt_sta_traces(this_traces,filt_param,varargin)
METHOD = 'midandmean';
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'METHOD')
        METHOD = varargin{v+1};
    end
end
% median and movmean along every row
this_filt_traces = this_traces;
for t = 1:size(this_traces,1)
    switch METHOD
        case 'midandmean'
            this_filt_traces(t,:) = movmean( medfilt1(this_traces(t,:),filt_param),filt_param);
        case 'lowess'
            % Bruno (Baker), 2020 used this filter, param = 0.002;
            % 0.01 worked better for our data.
            this_filt_traces(t,:) = smooth(this_traces(t,:),filt_param,'lowess');
    end
end
end

%% test
% figure;
% subplot(3,1,1)
% plot(this_traces')
% subplot(3,1,2)
% plot(this_filt_traces')
% subplot(3,1,3); hold on
% plot(this_filt_traces(1,:));plot(this_traces(1,:))