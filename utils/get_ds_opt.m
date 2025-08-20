function [ds_opt] = get_ds_opt(opt,ds_factor)
% devide all # frames related fields in opt by ds_factor
fds = fields(opt);
frame_fds = fds(contains(fds,'_frames'));
ds_opt = opt;
ds_opt.ds_factor = ds_factor;
for i = 1:numel(frame_fds)
    fd = frame_fds{i};
    ds_opt.(fd) = round(opt.(fd)/ds_factor);
end
ds_opt.frame_rate =  round(opt.frame_rate/ds_factor);
end

