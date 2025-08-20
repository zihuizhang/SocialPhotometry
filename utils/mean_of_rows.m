function [output] = mean_of_rows(traces,row_groups)
% average rows in traces that belong to the same 'row_group'
if iscell(row_groups)
    num_rows = numel(row_groups);
    output = nan(num_rows,size(traces,2));
    for r = 1:num_rows
        output(r,:) = nanmean(traces(row_groups{r},:),1);
    end
else
    num_rows = numel(unique(row_groups));
    output = nan(num_rows,size(traces,2));
    for r = 1:num_rows
        output(r,:) = nanmean(traces(row_groups==r,:),1);
    end

end
end