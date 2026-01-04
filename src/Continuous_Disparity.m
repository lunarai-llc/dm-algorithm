% written on 11/21/2018

% matlab tutorial for creating packages: 
% https://www.mathworks.com/help/matlab/matlab_prog/create-and-share-custom-matlab-toolboxes.html

% This is the discrete disparity method, specifically it implements the
% HD and NED method

% This includes the Binomial, Poisson, Geometric, Negative Binomial, 
% Poisson-Gamma, Poisson-Lognormal model

% The initial estimate is set to me moment estimate

% The optimization method is quasi newton

% Usage: Discrete_Disparity(X, 'HD', 'poiss')

% Optional argument is positional

% generate random variables
% normal:  X = normrnd(2, 1, 100, 1);
% Logistic: X = random('Logistic', 2, 1, 100, 1);
% t: X = trnd(nu, 100, 1); nu is the degree of freedom
% exponential: X = exprnd(.5, 100, 1);
% rayleigh: X = raylrnd(2, 100, 1);
% chisquare: X = chi2rnd(5, 100, 1);

function [ttheta] = Continuous_Disparity(Y, model, varargin)

defaultBand = [];
defaultDisparities = 'hd';
expectedDisparities = {'hd', 'ned'};
defaultKernels = 'epanechnikov';
expectedKernels = {'epanechnikov', 'normal', 'box', 'triangle'};
expectedModels = {'normal', 'gamma', 'beta', 'lognormal', 'weibull',...
                  'logistic' ,'t', 'exp', 'rayleigh', 'nakagami', ...
                  'rician', 'chisq', 'inversegaussian'};
              
inp = inputParser;
% validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarPosNum1 = @(x) isnumeric(x) ;
validScalarPosNum2 = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(inp, 'Y', validScalarPosNum1);
addRequired(inp, 'model', @(x) any(validatestring(x,expectedModels)));
addOptional(inp, 'band', defaultBand, validScalarPosNum2);
addParameter(inp, 'disparities', defaultDisparities, ...
                @(x) any(validatestring(x,expectedDisparities)) );
addParameter(inp, 'kernels', defaultKernels, ...
                @(x) any(validatestring(x,expectedKernels)) );             
% argin.addParamValue('geeMaxIter', 100, @(x) isnumeric(x) && x>0);            

parse(inp, Y, model, varargin{:});

Y = inp.Results.Y;
model = inp.Results.model;
band  = inp.Results.band;
disparities = inp.Results.disparities;
kernels = inp.Results.kernels;

if (isempty(band))
    if strcmp(kernels, 'normal')
        pd_kernel = fitdist(Y, 'Kernel','Kernel','normal');
    elseif strcmp(kernels, 'epanechnikov')
        pd_kernel = fitdist(Y, 'Kernel','Kernel','epanechnikov');
    elseif strcmp(kernels, 'box')
        pd_kernel = fitdist(Y, 'Kernel','Kernel','box');
    elseif strcmp(kernels, 'triangle')
        pd_kernel = fitdist(Y, 'Kernel','Kernel', 'triangle');
    end
else
%     cm = .55;
%     sm = 1.48*median(abs(X-median(X)));
%     band = cm*sm;
    pd_kernel = fitdist(Y, 'Kernel','Kernel','normal', 'BandWidth', band);        
end

% pd_kernel = fitdist(X, 'Kernel','Kernel','normal', 'BandWidth', h);

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');  
% options = optimoptions(@fminunc, 'Algorithm','quasi-newton');  

% disparity method
if strcmp(disparities,'hd')
    if strcmp(model, 'normal')
        obj_fun = @(uu, mu1, sigma1) (pdf(pd_kernel, uu).*normpdf(uu, mu1, sigma1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx(1), xx(2)), -Inf, Inf, 'ArrayValued', true);
        x0 = [mean(Y) std(Y)];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'gamma')
        obj_fun = @(uu, alpha1, beta1) (pdf(pd_kernel, uu).*gampdf(uu, alpha1, beta1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, Inf, 'ArrayValued', true);
        x0 = [mean(Y)^2/var(Y) var(Y)/mean(Y) ];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'beta')
        obj_fun = @(uu, alpha1, beta1) (pdf(pd_kernel, uu).*betapdf(uu, alpha1, beta1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, 1, 'ArrayValued', true);
        x0 = betafit(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'lognormal')
        obj_fun = @(uu, mu1, sigma1) (pdf(pd_kernel, uu).*lognpdf(uu, mu1, sigma1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, Inf, 'ArrayValued', true);
        x0 = lognfit(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'weibull')
        obj_fun = @(uu, a1, b1) (pdf(pd_kernel, uu).*wblpdf(uu, a1, b1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, Inf, 'ArrayValued', true);
        x0 = wblfit(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'logistic')
        obj_fun = @(uu, mu1, sigma1) (pdf(pd_kernel, uu).*pdf('Logistic', uu, mu1, sigma1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx(1), xx(2)), -Inf, Inf, 'ArrayValued', true);
        x0 = mle(Y, 'distribution', 'Logistic');
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 't')
        obj_fun = @(uu, nu1) (pdf(pd_kernel, uu).*tpdf(uu, nu1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), -Inf, Inf, 'ArrayValued', true);
        x0 = 2*var(Y)/(var(Y)-1);
        ttheta = fminunc(fun1, x0, options);    
        ttheta = round(ttheta);
    elseif strcmp(model, 'exp')
        obj_fun = @(uu, mu1) (pdf(pd_kernel, uu).*exppdf(uu, mu1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = 1/mean(Y);
        ttheta = fminunc(fun1, x0, options);     
    elseif strcmp(model, 'rayleigh') 
        obj_fun = @(uu, b1) (pdf(pd_kernel, uu).*raylpdf(uu, b1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = sqrt(.5*mean(Y.^2));
        ttheta = fminunc(fun1, x0, options);   
    elseif strcmp(model, 'nakagami')
        obj_fun = @(uu, mu1, omega1) (pdf(pd_kernel, uu).*pdf('Nakagami', uu, mu1, omega1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = mle(Y, 'Nakagami');
        ttheta = fminunc(fun1, x0, options);  
    elseif strcmp(model, 'rician')
        obj_fun = @(uu, s1, sigma1) (pdf(pd_kernel, uu).*pdf('Rician', uu, s1, sigma1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = mle(Y, 'Rician');
        ttheta = fminunc(fun1, x0, options);     
    elseif strcmp(model, 'chisq')
        obj_fun = @(uu, nu1) (pdf(pd_kernel, uu).*chi2pdf(uu, nu1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = mean(Y);
        ttheta = fminunc(fun1, x0, options);      
        ttheta = round(ttheta);
    elseif strcmp(model, 'inversegaussian')
        obj_fun = @(uu, mu1, lambda1) (pdf(pd_kernel, uu).* ...
                                pdf('InverseGaussian', uu, mu1, lambda1)).^.5;
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = mle(Y, 'InverseGaussian');
        ttheta = fminunc(fun1, x0, options);      
    else
        error('model not recogonized.');
    end
elseif strcmp(disparities,'ned')
    if strcmp(model, 'normal') % moment estimator for initial value
%         obj_fun = @(uu, mu1, sigma1) exp(-normpdf(uu, mu1, sigma1)./pdf(pd_kernel, uu)) ...
%                                      .*pdf(pd_kernel, uu);   
        obj_fun = @(uu, mu1, sigma1) exp(-pdf(pd_kernel, uu)./normpdf(uu, mu1, sigma1)) ...
                                         .*normpdf(uu, mu1, sigma1);   
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx(1), xx(2)), -15, 15, 'ArrayValued', true);
        x0 =  [mean(Y) std(Y)];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'gamma') % moment estimator for initial value
        obj_fun = @(uu, alpha1, beta1) exp(-pdf(pd_kernel, uu)./gampdf(uu, alpha1, beta1)) ...
                                           .*gampdf(uu, alpha1, beta1);   
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, 10, 'ArrayValued', true);
        x0 = [mean(Y)^2/var(Y) var(Y)/mean(Y)];
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'beta') % MLE for initial value
        obj_fun = @(uu, alpha1, beta1) exp(-pdf(pd_kernel, uu)./betapdf(uu, alpha1, beta1)) ...
                                           .*betapdf(uu, alpha1, beta1);   
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, 1, 'ArrayValued', true);
        x0 = betafit(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'lognormal') % moment estimator initial value
        obj_fun = @(uu, mu1, sigma1) exp(-pdf(pd_kernel, uu)./lognpdf(uu, mu1, sigma1)) ...
                                         .*lognpdf(uu, mu1, sigma1);   
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, Inf, 'ArrayValued', true);
        x0 = lognfit(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'weibull') % MLE for initial value
        obj_fun = @(uu, a1, b1) exp(-pdf(pd_kernel, uu)./wblpdf(uu, a1, b1)) ...
                                    .*wblpdf(uu, a1, b1);  
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx(1), xx(2)), 0, 20, 'ArrayValued', true);
        x0 = wblfit(Y);
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 'logistic') % MLE for initial value
        obj_fun = @(uu, mu1, sigma1) exp(-pdf(pd_kernel, uu)./pdf('Logistic', uu, mu1, sigma1)) ...
                                         .*pdf('Logistic', uu, mu1, sigma1);        
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx(1), xx(2)), -15, 15, 'ArrayValued', true);
        x0 = mle(Y, 'distribution', 'Logistic');
        ttheta = fminunc(fun1, x0, options);
    elseif strcmp(model, 't') % moment
        obj_fun = @(uu, nu1) exp(-pdf(pd_kernel, uu)./tpdf(uu, nu1)) ...
                                 .*tpdf(uu, nu1);      
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx), -15, 15, 'ArrayValued', true);
        x0 = 2*var(Y)/(var(Y)-1);
        ttheta = fminunc(fun1, x0, options); 
        ttheta = round(ttheta);
    elseif strcmp(model, 'exp') % moment
        obj_fun = @(uu, mu1) exp(-pdf(pd_kernel, uu)./exppdf(uu, mu1)) ...
                                 .*exppdf(uu, mu1);   
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx), 0, 20, 'ArrayValued', true);
        x0 = 1/mean(Y);
        ttheta = fminunc(fun1, x0, options);     
    elseif strcmp(model, 'rayleigh') % moment
        obj_fun = @(uu, b1) exp(-pdf(pd_kernel, uu)./raylpdf(uu, b1)) ...
                                 .*raylpdf(uu, b1);   
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx), 0, 20, 'ArrayValued', true);
        x0 = sqrt(.5*mean(Y.^2));
        ttheta = fminunc(fun1, x0, options);       
    elseif strcmp(model, 'nakagami') % mle
        obj_fun = @(uu, mu1, omega1) exp(-pdf(pd_kernel, uu)./pdf('Nakagami', uu, mu1, omega1)) ...
                                         .*pdf('Nakagami', uu, mu1, omega1);           
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx), 0, 20, 'ArrayValued', true);
        x0 = mle(Y, 'Nakagami');
        ttheta = fminunc(fun1, x0, options);        
    elseif strcmp(model, 'rician') % mle
        obj_fun = @(uu, s1, sigma1) exp(-pdf(pd_kernel, uu)./pdf('Rician', uu, s1, sigma1)) ...
                                         .*pdf('Rician', uu, s1, sigma1);  
        fun1 = @(xx) -integral(@(uu)obj_fun(uu, xx), 0, 20, 'ArrayValued', true);
        x0 = mle(Y, 'Rician');
        ttheta = fminunc(fun1, x0, options);     
    elseif strcmp(model, 'chisq')
        obj_fun = @(uu, nu1) exp(-pdf(pd_kernel, uu)./chi2pdf(uu, nu1)) ...
                                  .*chi2pdf(uu, nu1);  
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx), 0, 20, 'ArrayValued', true);
        x0 = mean(Y);
        ttheta = fminunc(fun1, x0, options);      
        ttheta = round(ttheta);  
    elseif strcmp(model, 'inversegaussian')
        obj_fun = @(uu, mu1, lambda1) exp(-pdf(pd_kernel, uu)./ ...
                                        pdf('InverseGaussian', uu, mu1, lambda1)) ...
                                        .*pdf('InverseGaussian', uu, mu1, lambda1);          
        fun1 = @(xx) integral(@(uu)obj_fun(uu, xx), 0, Inf, 'ArrayValued', true);
        x0 = mle(Y, 'InverseGaussian');
        ttheta = fminunc(fun1, x0, options);          
    else
        error('model not recogonized.');
    end
else
    error('disparity not recogonized.');
end

end



% argin.addRequired('Y', @isnumeric);
% argin.addRequired('disparity', @(x) strcmpi(x,'hd') || ...
%     strcmpi(x,'ned'));
% argin.addRequired('model', @(x) strcmpi(x,'normal') || ...
%     strcmpi(x,'gamma') || strcmpi(x,'beta') || strcmpi(x,'lognormal') || ...
%     strcmpi(x,'weibull') || strcmpi(x,'logistic') || strcmpi(x,'t') || ...
%     strcmpi(x, 'F') || strcmpi(x, 'chisq') || strcmpi(x, 'exp') || ...
%     strcmpi(x, 'dexp') );
% argin.addRequired('ker', @(x) strcmpi(x,'normal') || ...
%     strcmpi(x,'epanechnikov') || strcmpi(x, 'box') || strcmpi(x, 'triangle'));
% argin.addParameter('band', [], @numeric);

% argin.addParamValue('maxiter', 1000, @(x) isnumeric(x) && x>0);

% argin.addParameter('aa', 1, @(x) isnumeric(x) && x>0);

% this is critical: to pass optional parameter into use
% argin.parse(Y, disparity, varargin{:});
% aa = argin.Results.aa;
% 
% % set default value
% disparity = upper(disparity);
% if (strcmp(disparity, 'HD'))
%     if (isempty(aa))
%         aa = 1;       
%     end
% end

% sample size
% n = length(Y);
