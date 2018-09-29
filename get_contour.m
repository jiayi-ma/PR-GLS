function contour = get_contour(mask, K)

[y,x] = ind2sub(size(mask), find(mask));
contour = bwtraceboundary(mask, [y(1),x(1)], 'N', 8, inf, 'clockwise');
contour = contour(:,end:-1:1);
if nargin == 2
    n = size(contour,1);
    g = linspace(1, n, K+1);
    p = repmat(g' - floor(g'), [1,2]);
    contour = contour(floor(g),:).*(1-p) + contour(mod(floor(g),n)+1,:).*p;
    contour(end,:) = [];
end