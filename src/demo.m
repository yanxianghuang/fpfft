% demo on compare the errors in (2048-point) fixed-point FFT
% created by Yanxiang on 16, Spet, 2016
% yanxiang.huang@imec.be

clear; close all;
np = 2048; % 2k fft test
systematic_comparison = true;


%% 1. input
figure(1);
in = ((rand(1,np)-0.5).*2);
%in = (sin(1:np));
subplot(5,1,1);
plot(in);
title('input');


input_frac = 11; %% It is tested the result is correct until input_frac = 28, which is way sufficent
tw_frac = input_frac-5; %% It is tested the result is correct until tw_frac = 28, which is way sufficient










%% 2. my fixed_point, using the generated c-mex file

% IMPORTANT SIZE INFORMATION !!!!!!!!!!!!!!!!!!
% the bit_width of input is calculated as input_frac + input_integer
% i) input_integer is :  " ceil(  log2( max(abs(real(in))) +1e-3  )  )  +1 "
% so for QAM_256, the maximum value is 15, therefore input_integer is
% log2(15) + 1 = 5
% ii) in the example above, input are in [1, -1], so input_integer = 2;
% so the realistic-size of you data are input_frac + 2;

% similiar is for the twidder factor of FFT, However, as the range of tw is
% always [-1, 1], the tw_integer is 2 ( log2(1+1e-3) + 1; note the +1e-3),
% so its size is always tw_frac+2

tic;
out1 = fix_fft2k(in, input_frac, tw_frac); % easy to use: complex input, input_frac_size, tw_frac size
toc;
subplot(5,1,2);
plot( fftshift(abs(fft(in)) ));
title('Perfect FFT output');

subplot(5,1,3);
plot( fftshift(abs(out1) ));
title('My FFT output');

subplot(5,1,4);
my_error = fftshift( abs(out1)-abs(fft(in)) );
plot( my_error );
title('My FFT error');
fprintf(sprintf('###My C-mex FFT rms error %f \n', rms(my_error)));





%% 3. error analysis (optional)
if systematic_comparison,
  my_error_grid = 0;
  for the_input_frac = 1:1:30,
    for the_tw_frac = 1:1:26,
      out1 = fix_fft2k(in, the_input_frac, the_tw_frac);
      my_error_grid(the_input_frac, the_tw_frac) = rms( fftshift( abs(out1)-abs(fft(in)) ) );
    end
  end
  figure(2);
  contour( log10(my_error_grid) ); %% in log10 domain
  title('my fixed-point FFT error vs. input-size and tw-size');
  xlabel('tw-frac size');
  ylabel('input-frac size');
  fprintf('You can find that in figure(2), with lines showing same error-level, the corners are almost always at input_frac = tw_frac+5. \nSo keep that trend as the rule-of-thumb in your design.\n');
end  
  
  
  
  
%% 4. quantize only inputs, optimistic
tic;
fixed_in = sfi(in, input_frac+2);
fixed_in = fixed_in.data;

out2 = fft(fixed_in);
toc;
figure(1),subplot(5,1,5);
fi_error = fftshift( abs(out2)-abs(fft(in)) );
plot( fi_error );
title('the only-quantize input FFT error');
fprintf(sprintf('###Only-quantize-input built-in FFT rms error %f \n', rms(fi_error)));

fi_error_grid = 0;
for the_input_frac = 1:1:30,
  fixed_in = sfi(in, the_input_frac+2);
  fixed_in = fixed_in.data;
  out2 = fft(fixed_in);
  fi_error_grid(the_input_frac) = rms( fftshift( abs(out2)-abs(fft(in)) ) );
end

figure(3);
semilogy( fi_error_grid);
title('if only quantization input, FFT error decrease as input-size increase');
xlabel('frac_width');ylabel('rms error');
