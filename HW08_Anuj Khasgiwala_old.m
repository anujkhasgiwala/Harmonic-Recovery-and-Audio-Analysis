%%% =======================================================================
%%% HW08: Harmonic Recovery and Audio Analysis
%%% Author: Anuj Khasgiwala
%%% =======================================================================
function problem_01()
  warning("off", "Octave:broadcast");
  DELTA = (2*pi)/6283;
  threshold = 0.1;
  t = -pi:DELTA:pi;

  %% combined sinusoid
  data = load("D:/workspace/USU-Assignments/CS 6810 - Wavelets and Wavelet Algorithms/Assignment 8/hw08/data.txt");
  data = data'(:)';

  
  f = 1:150;
  sine_coeffs = [];
  cos_coeffs = [];
  for fi=f
      sine_coeff = 1/pi*sum(data.*sin(fi*t))*DELTA;
      if ( abs(sine_coeff) > threshold )
        disp(strcat(strcat(strcat("a",num2str(fi))," = "), num2str(abs(sine_coeff))));
        sine_coeffs(fi) = abs(sine_coeff);
      else
        sine_coeffs(fi) = 0.0;
      endif
  %% ================= Recovering cosine Coefficients ======================
      cos_coeff = 1/pi*sum(data.*cos(fi*t))*DELTA;
      if ( abs(cos_coeff) > threshold )
        disp(strcat(strcat(strcat("b",num2str(fi))," = "), num2str(abs(cos_coeff))));
        cos_coeffs(fi) = abs(cos_coeff);
      else
        cos_coeffs(fi) = 0.0;
      endif
  end

  %% plot sine coefficients
  figure;
  plot(1:150, sine_coeffs);
  xlabel('Freq');
  ylabel('Sine Coeff');
  title('Present Sine Coefficients');
endfunction