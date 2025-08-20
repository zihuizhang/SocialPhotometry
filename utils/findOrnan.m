function [idx,ERROR] = findOrnan(a,b)
idx = find(a>b,1);
ERROR = 0;
if isempty(idx)
    idx = 1;
    ERROR = 1;
end

end

