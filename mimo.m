clear; clc; close all;

% Parameters
numTx = 4;    % Number of transmitters 
numRx = 4;    % Number of receivers 
fc = 100e3;   % Carrier frequency of the sensor
c = 1500;     % Speed of sound in water 
lambda = c / fc;  % Wavelength of the ultrasonic 
d = lambda / 2;   % Spacing between MIMO elements 
N = 1000;     % Number of time samples
fs = 1e6;     % Sampling frequency (1 MHz)

% Generate the transmitted signal (sinusoidal wave)
t = (0:N-1) / fs;   % Time 
signal = cos(2 * pi * fc * t); % Transmit signal , can also use sine for other results

% Create MIMO transmit and receive arrays
tx_signal = repmat(signal, numTx, 1);  % Repeat the signal across all transmitters
rx_signal = zeros(numRx, N);           % received signal matrix

% Simulate signal propagation and noise at each receiver
for rxIdx = 1:numRx
    rx_signal(rxIdx, :) = tx_signal(1, :) + 0.1 * randn(1, N);  % Add noise to each received signal
end

% % Create a MIMO channel model with beamforming
% tx_array = phased.ULA('NumElements',numTx,'ElementSpacing',d); % Transmitter array
% rx_array = phased.ULA('NumElements',numRx,'ElementSpacing',d); % Receiver array
% 
% % Apply beamforming towards the target
% beamformer = phased.PhaseShiftBeamformer('SensorArray',rx_array, ...
%     'OperatingFrequency',fc,'Direction',[theta_target; 0], 'WeightsOutputPort', true);
% 
% % Transmit the signal from each element (MIMO)
% tx_signal = repmat(signal, numTx, 1); % Transmit signal on each antenna
% 
% % Simulate the received signal with some noise and interference
% interference_signal = 0.5 * cos(2 * pi * fc * t + pi/4); % Interference from another source
% noise = 0.1 * randn(size(tx_signal)); % Additive Gaussian noise
% 
% % Set the angles for the received signal (must be repeated for each sample)
% theta_target_vec = repmat([theta_target; 0], 1, N);
% theta_interference_vec = repmat([theta_interference; 0], 1, N);


% Visualization
figure;
subplot(2,1,1);
plot(t*1e3, tx_signal(1,:)); %tx + noise
title('Transmitted Signal');
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;

subplot(2,1,2);
plot(t*1e3, rx_signal(1,:)); %rx + nnoise
title('Received Signal (with Noise)');
xlabel('Time (ms)');
ylabel('Amplitude');
grid on;

% SNR
snr_value = snr(rx_signal(1,:), 0.1 * randn(1, N));  % SNR for one received signal
fprintf('SNR: %.2f dB\n', snr_value);
