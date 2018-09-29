function [P, P1, Pt1, PX, L]=compute_P(X,T, sigma2 ,outlier, s_matrix)

[N, D] = size(X);
M = size(T, 1);
% s_matrix = ones(N, M) / M;
if D == 2
    a = (max(X(:,1)) - min(X(:,1)))*(max(X(:,2)) - min(X(:,2)));
    b = (max(T(:,1)) - min(T(:,1)))*(max(T(:,2)) - min(T(:,2)));
else
    a = (max(X(:,1)) - min(X(:,1)))*(max(X(:,2)) - min(X(:,2)))*(max(X(:,3)) - min(X(:,3)));
    b = (max(T(:,1)) - min(T(:,1)))*(max(T(:,2)) - min(T(:,2)))*(max(T(:,3)) - min(T(:,3)));
end
a = max(a, b);
ksig = -2.0 * sigma2;
outlier_tmp=outlier*(-ksig*3.14159265358979)^(0.5*D)/((1-outlier)*a);

P = repmat(X,[1 1 M])-permute(repmat(T,[1 1 N]),[3 2 1]);
P = squeeze(sum(P.^2,2));
P = P/ksig;
P = exp(P).*s_matrix;
sum_tmp = sum(P,2) + outlier_tmp;
P = P./repmat(sum_tmp, 1, M);
P = P';
P1 = P*ones(N,1);
Pt1 = P'*ones(M,1);
PX = P*X;
L = sum(log(sum_tmp)) + D*N*log(sigma2)/2;