% Problem 1
dataFile = fopen('data.txt','r');
sinusoid_data = fscanf(dataFile,'%f',[1 6284]);
fclose(dataFile);

% Recovering Sine Coefficients
trig = 0;
f = 1:150;
DELTA = 0.001;
t = -pi:DELTA:pi;
coeffs = zeros(size(f));
for fi = f
    if(trig == 0)
        coeff = 1/pi*sum(sinusoid_data.*sin(fi*t))*DELTA;
    else
        coeff = 1/pi*sum(sinusoid_data.*cos(fi*t))*DELTA;
    end
    if(abs(coeff) > 0.001)
        coeffs(fi) = coeff;
        if (abs(coeff) > 0.03)
            if(trig == 0)
                display(['b',num2str(fi), ' = ', num2str(coeff)]);
            else
                display(['a',num2str(fi), ' = ', num2str(coeff)]);
            end
        end
    end
end
fprintf('\n');
% plot
figure;
plot(f, coeffs);
xlabel('Freq');
ylabel('Coeff');

% Recovering Cosine Coefficiens
trig = 1;
f = 1:150;
DELTA = 0.001;
t = -pi:DELTA:pi;
coeffs = zeros(size(f));
for fi = f
    if(trig == 0)
        coeff = 1/pi*sum(sinusoid_data.*sin(fi*t))*DELTA;
    else
        coeff = 1/pi*sum(sinusoid_data.*cos(fi*t))*DELTA;
    end
    if(abs(coeff) > 0.001)
        coeffs(fi) = coeff;
        if (abs(coeff) > 0.03)
            if(trig == 0)
                display(['b',num2str(fi), ' = ', num2str(coeff)]);
            else
                display(['a',num2str(fi), ' = ', num2str(coeff)]);
            end
        end
    end
end
fprintf('\n');
% plot
figure;
plot(f, coeffs);
xlabel('Freq');
ylabel('Coeff');

% Problem 2
% file1
[y,Fs] = audioread(strcat('ODE_TO_JOY_BOTH_HANDS','.wav'));
Nsamps = length(y);

y_fft = abs(fft(y));
y_fft = y_fft(1:Nsamps/2);
f = Fs*(0:Nsamps/2-1)/Nsamps;

figure;
plot(f, y_fft);
xlim([0 1000]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['Frequency Response of ',strrep('ODE_TO_JOY_BOTH_HANDS','_',' ')]);
disp(['Top 5 Frequencies in ','ODE_TO_JOY_BOTH_HANDS']);
disp('1 - E4');
disp('2 - G4');
disp('3 - D4');
disp('4 - C4');
disp('5 - D5');
fprintf('\n');

% file2
[y,Fs] = audioread(strcat('ODE_TO_JOY_ORCHESTRA','.wav'));
Nsamps = length(y);

y_fft = abs(fft(y));
y_fft = y_fft(1:Nsamps/2);
f = Fs*(0:Nsamps/2-1)/Nsamps;

figure;
plot(f, y_fft);
xlim([0 1000]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['Frequency Response of ',strrep('ODE_TO_JOY_ORCHESTRA','_',' ')]);
disp(['Top 5 Frequencies in ','ODE_TO_JOY_BOTH_HANDS']);
disp('1 - D5');
disp('2 - G5');
disp('3 - A5');
disp('4 - E5');
disp('5 - A4');
fprintf('\n');

% file3
[y,Fs] = audioread(strcat('ODE_TO_JOY_RIGHT_HAND','.wav'));
Nsamps = length(y);

y_fft = abs(fft(y));
y_fft = y_fft(1:Nsamps/2);
f = Fs*(0:Nsamps/2-1)/Nsamps;

figure;
plot(f, y_fft);
xlim([0 1000]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['Frequency Response of ',strrep('ODE_TO_JOY_RIGHT_HAND','_',' ')]);
disp(['Top 5 Frequencies in ','ODE_TO_JOY_RIGHT_HAND']);
disp('1 - E4');
disp('2 - G4');
disp('3 - D4');
disp('4 - D5');
disp('5 - C5');
