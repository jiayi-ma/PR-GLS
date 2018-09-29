function  [P, C, W, iter, T] = prgls_GRBF(X, Y, beta, lambda, max_it, tol, viz, outliers, corresp,sigma2, t, nsc);

[N, D]=size(X); [M, D]=size(Y);

% Initialization
T=Y; iter=0;  ntol=tol+10; W=zeros(M,D);
if ~exist('sigma2','var') || isempty(sigma2) || (sigma2==0), 
    sigma2=(M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
end
sigma2_init=sigma2;

% Construct affinity matrix G
G=cpd_G(Y,Y,beta);

% configuration of SC
mean_dist_global=[]; % use [] to estimate scale from the data
nbins_theta=12;
nbins_r=5;
nsamp1=size(X,1);
nsamp2=size(Y,1);
ndum1=max(0, nsamp2-nsamp1);
ndum2=max(0, nsamp1-nsamp2);
eps_dum=0.15;
r_inner=1/8;
r_outer=2;
n_iter=30;
out_vec_1=zeros(1,nsamp1); 
out_vec_2=zeros(1,nsamp2);

% c1=mean(X);
% t1=atan2(X(:,2)-c1(2), X(:,1)-c1(1));
% [BH1,mean_dist_1]=sc_compute(X',t1',mean_dist_global,nbins_theta,...
%     nbins_r,r_inner,r_outer,out_vec_1);
[BH1,mean_dist_1]=sc_compute(X',zeros(1,nsamp1),mean_dist_global,nbins_theta,...
    nbins_r,r_inner,r_outer,out_vec_1);

iter=0; ntol=tol+10; L=1;
while (iter<max_it) && (ntol > tol) && (sigma2 > 1e-8) %(sigma2 > 1e-8)
%while (iter<max_it)  && (sigma2 > 1e-8) %(sigma2 > 1e-8)
    
    if mod(iter, nsc)==0
%     c2=mean(T);
%     t2=atan2(T(:,2)-c2(2), T(:,1)-c2(1));
%     [BH2,mean_dist_2]=sc_compute(T',t2',mean_dist_1,nbins_theta,...
%         nbins_r,r_inner,r_outer,out_vec_2);
    [BH2,mean_dist_2]=sc_compute(T',zeros(1,nsamp2),mean_dist_1,nbins_theta,...
        nbins_r,r_inner,r_outer,out_vec_2);
    % compute pairwise cost between all shape contexts
    costmat=hist_cost_2(BH1,BH2);
    % pad the cost matrix with costs for dummies
    nptsd=nsamp1+ndum1;
    costmat2=eps_dum*ones(nptsd,nptsd);
    costmat2(1:nsamp1,1:nsamp2)=costmat;
%     disp('running hungarian alg.')
    cvec=hungarian(costmat2);
%     cvec = 1:98;
%     disp('done.')
    
    P_prior = zeros(N, M);
    
    if N > M
        for i = 1:M
            P_prior(cvec(i),i) = t;
        end
    else
        for i = 1:N
            P_prior(i,cvec(i)) = t;
        end
    end
    P_prior = P_prior + (1-t)/M;
    idx = sum(P_prior,2) < 1;
    P_prior(idx,:) = 1/M;
    end
    
    L_old=L;
    
    [P, P1, Pt1, PX, L]=compute_P(X,T, sigma2 ,outliers, P_prior); st='';
    
    disp([' PR-GLS nonrigid ' st ' : dL= ' num2str(ntol) ', iter= ' num2str(iter) ' sigma2= ' num2str(sigma2)]);
    
    L=L+lambda/2*trace(W'*G*W);
    ntol=abs((L-L_old)/L);

    dP=spdiags(P1,0,M,M); % precompute diag(P)
    W=(dP*G+lambda*sigma2*eye(M))\(PX-dP*Y);

    % update Y postions
    T=Y+G*W;

    Np=sum(P1);sigma2save=sigma2;
    sigma2=abs((sum(sum(X.^2.*repmat(Pt1,1,D)))+sum(sum(T.^2.*repmat(P1,1,D))) -2*trace(PX'*T)) /(Np*D));
    
    outliers = 1 - Np/N; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if outliers>0.99
        outliers = 0.99;
    elseif outliers<0.01
        outliers = 0.01;
    end
    
%     outliers = 0.01;

    % Plot the result on current iteration
    if viz, figure(1); cpd_plot_iter(X, T); end;
    iter=iter+1;

end


disp('PR-GLS registration succesfully completed.');

%Find the correspondence, such that Y corresponds to X(C,:)
if corresp, C=cpd_Pcorrespondence(X,T,sigma2save,outliers); else C=0; end;
