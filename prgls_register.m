function [Transform, C]=prgls_register(X, Y, opt)

[M,D]=size(Y); [N, D2]=size(X);
% Check the input options and set the defaults
if nargin<3, error('Not enough input parameters.'); end;
if ~isfield(opt,'normalize') || isempty(opt.normalize), opt.normalize = 1; end;
if ~isfield(opt,'max_it') || isempty(opt.max_it), opt.max_it = 100; end;
if ~isfield(opt,'tol') || isempty(opt.tol), opt.tol = 1e-5; end;
if ~isfield(opt,'viz') || isempty(opt.viz), opt.viz = 1; end;
if ~isfield(opt,'corresp') || isempty(opt.corresp), opt.corresp = 0; end;
if ~isfield(opt,'outliers') || isempty(opt.outliers), opt.outliers = 0.1; end;
if ~isfield(opt,'sigma2') || isempty(opt.sigma2), opt.sigma2 = 0; end;
if ~isfield(opt,'sparse') || isempty(opt.sparse), opt.sparse = 0; end;
if ~isfield(opt,'beta') || isempty(opt.beta), opt.beta = 2; end;
if ~isfield(opt,'lambda') || isempty(opt.lambda), opt.lambda = 3; end;
if ~isfield(opt,'nsc') || isempty(opt.nsc), opt.nsc = 10; end;

% checking for the possible errors
if D~=D2, error('The dimension of point-sets is not the same.'); end;
if (D>M)||(D>N), disp('The dimensionality is larger than the number of points. Possibly the wrong orientation of X and Y.'); end;
if (D<=1) || (D>3), opt.viz=0; end;

% check if mex functions are compiled yet
if ~exist('cpd_P','file')
    disp('Looks like you have not compiled CPD mex files yet (needs to be done once)');
    disp('Running cpd_make.m for you ...'); tic;
    cpd_make;
end

% Convert to double type, save Y
X=double(X);  
Y=double(Y); Yorig=Y; 

% default mean and scaling
normal.xd=0; normal.yd=0;
normal.xscale=1; normal.yscale=1;

% Normalize to zero mean and unit variance
if opt.normalize, [X,Y,normal]=cpd_normalize(X,Y); end;

% disp(['%%%%% Starting CPD-' upper(opt.method) ' registration. %%%' ]); tic;

if opt.sparse
    [P, C, W, iter, T] = prgls_GRBF_sparse(X, Y, opt.beta, opt.lambda, opt.max_it, opt.tol, opt.viz, opt.outliers, opt.corresp, opt.sigma2, opt.t, opt.nsc); 
else
    [P, C, W, iter, T] = prgls_GRBF(X, Y, opt.beta, opt.lambda, opt.max_it, opt.tol, opt.viz, opt.outliers, opt.corresp, opt.sigma2, opt.t, opt.nsc); 
end

Transform.iter=iter;
Transform.Y=T;
Transform.normal=normal;
Transform.P=P;
% Transform.W = W(:);

Transform.beta=opt.beta;
Transform.W=W;
Transform.Yorig=Yorig;
Transform.s=1;
Transform.t=zeros(D,1);
if opt.normalize,
     Transform.Y=T*normal.xscale+repmat(normal.xd,M,1);
end 






