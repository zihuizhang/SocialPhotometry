function [ peak_frame ] = template_matching( trace,template,Th )
% implementing scaled template matching algorithm
% J. Clements Detection of spontaneous synaptic events with an optimally scaled
% template 1997

N = length(template);
sum_template = sum(template);
sq_template = sum(template.^2);
scale_denom =  sq_template - sum_template*sum_template/N;

for i = 1:length(trace)- N
    this_data = trace(i:i+N-1);
    sum_data = sum(this_data);
    scale(i) = (sum(template.*this_data)-sum_template*sum_data/N)/scale_denom;
    offset(i) = (sum_data-scale(i)*sum_template)/N;
    
    fitted_temp = template.*  scale(i) +offset(i);
    this_sse = sum((this_data - fitted_temp).^2);
    std_error = sqrt(this_sse/(N-1));
    inter_detect_criterion(i)= scale(i)/std_error;
    
end

peak_frame = find(inter_detect_criterion>Th);





end

