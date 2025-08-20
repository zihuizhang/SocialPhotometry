function [robust_db_green,this_db_green] = debleach(this_green,opt)
func = @(x)FP_DEBLEACHED(x(opt.discard_first_frames+1:end)',120, opt.frame_rate, 100);
[this_db_green,robust_db_green] = func(this_green); %Get debleached trace; takes several sec 
% pad nans to start
this_db_green = [ nan(opt.discard_first_frames,1); double(this_db_green)];
robust_db_green = [ nan(opt.discard_first_frames,1); double(robust_db_green)];

end

