function [S_f,f,amp,pha,power,psd] = fftspectrum(s_t, nfft, fs, freqopt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [S_f,f,amp,pha,power,psd] = fftspectrum(s_t, nfft, fs, freqopt)
%
% This function would compute the discrete fast fourier transform, and return
% the one-sided or two-sided signal in freq domain, and corresponding freq 
% axis, amplitude spectrum, phase spectrum, power spectrum and power spectral
% density.
%
% See also 'testPSD.mlx', 'testfft.m'.
%
%
% Chao Song, chaosong@princeton.edu
% First created date:   2022/08/04
% Last modified date:   2022/08/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defval('freqrange','onesided');

S_f = fft(s_t, nfft);

if isequal(freqopt, 'onesided')
  S_f = S_f(1:nfft/2+1, :);  %discard half of points, needs to be calibrated later
  f = (0:nfft/2)*(fs/nfft); %frequency axis
  amp = abs(S_f)/nfft;   %magnitude spectrum
  amp(2:end-1, :) = 2*amp(2:end-1, :);  %times 2 since preserving only the positive freqs
  
  pha = angle(S_f);   %phase spectrum
  
  power = abs(S_f/nfft).^2;
  power(2:end-1, :) = 2*power(2:end-1, :);  %power spectrum == amp*2 /2, 'help periodogram'
  
  psd = abs(S_f).^2/nfft/fs;
  psd(2:end-1, :) = 2*psd(2:end-1, :);  %power spectral density
  
elseif isequal(freqopt, 'twosided')
  f = (0:nfft-1)*(fs/nfft); %frequency axis
  amp = abs(S_f)/nfft;   %magnitude spectrum
  pha = angle(S_f);
  power = abs(S_f/nfft).^2;
  psd = abs(S_f).^2/nfft/fs;
end
 


