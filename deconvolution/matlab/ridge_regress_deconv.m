function neuro = ridge_regress_deconv(BOLD,TR,alpha,NT)

% Deconvolves a preprocessed BOLD signal into neuronal time series 
% based on discrete cosine set and ridge regression.
% This function does not use confound regressors (e.g. motion).
% This function does not perform whitening and temporal filtering.
% The BOLD input signal must already be pre-processed.
%
% FORMAT neuro = ridge_regress_deconv(BOLD,TR,alpha,NT)
% BOLD  - Preprocessed BOLD signal
% TR    - Time repetition, [s]
% alpha - Regularization parameter (default: 0.005)
% NT    - Microtime resolution (number of time bins per scan, default: 16)
%
% _________________________________________________________________________
% Copyright (C) 2024 Ruslan Masharipov
% Contact email: masharipov@ihb.spb.ru

% Setup variables
if nargin < 2
    error('Define time repetition (TR)')
elseif nargin < 3
    alpha = 0.005;
    NT = 16;
elseif  nargin == 3
    NT = 16;
end

dt = TR/NT;                      % Length of time bin, [s]
N = length(BOLD);                % Scan duration, [dynamics] 
k = 1:NT:N*NT;                   % Microtime to scan time indices

% Create canonical HRF in microtime resolution (identical to SPM cHRF)
t = 0:dt:32;
hrf = gampdf(t,6) - gampdf(t,NT)/6;
hrf = hrf'/sum(hrf);

% Create convolved discrete cosine set
M = N*NT + 128;
n = (0:(M -1))';
xb = zeros(size(n,1),N);
xb(:,1) = ones(size(n,1),1)/sqrt(M);
for j=2:N
    xb(:,j) = sqrt(2/M)*cos(pi*(2*n+1)*(j-1)/(2*M));
end

Hxb = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:);

% Perform ridge regression
C = (Hxb'*Hxb + alpha*eye(length(Hxb)))\(Hxb'*BOLD);

% Recover neuronal signal
neuro = xb*C;

end