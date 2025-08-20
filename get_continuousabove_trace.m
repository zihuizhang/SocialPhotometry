function [output_trace,episodes] = get_continuousabove_trace(input_trace,mintime,ignorepts)
try
episodes = pt_continuousabove(input_trace,0,0.1,mintime,Inf,ignorepts);
output_trace = zeros(size(input_trace));
for t = 1:size(episodes,1)
    output_trace(episodes(t,1):episodes(t,2)) = 1;
end
catch
    output_trace = input_trace;
    episodes = [];
end
end