function y = nansem(y)

y = nanstd(y)./sqrt(sum(~isnan(y)));
y(y == 0) = nan;

end

