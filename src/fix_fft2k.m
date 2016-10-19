function [ output_args ] = fix_fft2k( x, x_frac, tw_frac )
%FIX_FFT2K Summary of this function goes here
%   Detailed explanation goes here
% input to FFT is x, x_frac and tw_frac meaning the number of bit to
% represent after the decimal_point.

% e.g. if x = 011.000001, then x_frac is 6;
% for tw, as they are between [-1, 1], then the integer_bit is always 2,
% and tw-bits are flexible.


NPoint = 2048;
NStage = 11; % log2(NPoint)


tw = exp( -2*pi*1i*(0:NPoint/2)/NPoint );

tw_2k_real = real(tw(1:end-1));
tw_2k_imag = imag(tw(1:end-1));



[real_out, imag_out ] = fix_fft(real(x), imag(x), tw_2k_real, tw_2k_imag,...
    int32(NPoint), int32(NStage), int32(x_frac), int32(tw_frac) );


%% we do not have bit reorder in FFT.cpp file, so add one here in
output_args = bitrevorder(real_out + 1j.*imag_out); 

end

