
clear; 
close; 
clc;

if 1
    initialization;
end

load save_fish_def_5_1.mat
X = x1;
Y = y2a;

opt.outliers = 0.5;
opt.viz = 1;
opt.t = 0.9;
opt.sparse = 0;
opt.nsc = 5;
% opt.normalize = 0;
% opt.beta = 2;
% opt.lambda = 3;
% opt.tol = 1e-10;

[Transform, C]=prgls_register(Y, X, opt);
V = Transform.Y;
       
figure,cpd_plot_iter(X, Y); axis off; title('Before');
figure,cpd_plot_iter(Transform.Y, y2a); axis off; title('After registering Y to X');
