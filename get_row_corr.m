function [shifts,correlation,cshifts,xcenter] = get_row_corr(A,B)
% get phase shift between each row in A and B
% A B must of same size
shifts = nan(1,size(A,1));
cshifts = shifts;
correlation = nan(size(A));
for i = 1:size(A,1)
    [c,lags] = xcorr(A(i,:),B(i,:),'normalized');
    if i == 1
        correlation = nan(size(A,1),length(c));
    end
    [~,Lag]=max(c);
    shifts(i) = lags(Lag); cshifts(i)= Lag;
    correlation(i,:) = c;
    xcenter = find(lags==0);
end

end

