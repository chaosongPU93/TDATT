function [S_f, S_f_conj] = time2freq(s_t, nfft)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [S_f, S_f_conj] = time2freq(s_t, nfft)
%
% This function would transform the time-domain signal to frequency domain,
% and also compute its complex conjugate. 
%
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/01/13
% Last modified date:   2022/01/13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the columns of signal
ncol = size(s_t,2);

%do for each column
S_f = zeros(nfft,ncol);
for j = 1: ncol
% fast fourier transform
  S_f(:, j) = fft(s_t(:, j), nfft);
  
end

% % fast fourier transform
% S_f = fft(s_t, nfft);

% complex conjugate
S_f_conj = conj(S_f);