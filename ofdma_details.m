clearvars;close all;clc;
SER =0;
% 16QAM MODULATION
M=16;
k = log2(M);
c= [-3 -1 1 3];
[a b] = meshgrid(c);
s = [reshape(a,[1,M]);reshape(b,[1,M])];
s = s(1,:)+sqrt(-1)*s(2,:);
K = 64; % number of ofdm carriers
CP = floor(K/4);%Length of the cyclic prefix: 25% of the block
P = 8;%number of pilor carriers per OFDM block
pilotValue = 3+3i; %the known value of each pilot transmits
%CALCULATING PILOT CARRIERS AND DATA CARRIERS
allcarriers = zeros(1,64);
for i =1:64
    allcarriers(i) = i;
end
%Adding pilot value at every 8th sample
pilotcarriers = allcarriers(1:floor(K/P):end);
%Let the last sample be pilot
pilotcarriers = [pilotcarriers allcarriers(end)];
%Remaining indices are deddicated to datacarriers 
datacarriers = setdiff(allcarriers,pilotcarriers);
%Visualizing Data Carriers and Pilot Carriers
figure('Renderer', 'painters', 'Position', [10 10 900 300])
plottin = zeros(1,64);
plot(pilotcarriers, plottin(pilotcarriers), 'b.','Markersize',20);
hold on;
plot(datacarriers, plottin(datacarriers), 'r.','Markersize',20);
legend('Pilot Carriers','Data Carriers')
xlabel('Carrier Index');
title('Location of Data Carriers and Pilot Carriers')
grid minor;
%16 QAM therefore 4 bits are required to represent each symbol
mu = 4;%BITS PER SYMBOL
payloadBits_per_OFDM = length(datacarriers)*mu; %number of payload bits per OFDM symbol
%GENERATING BIT STREAM
Txbits = randi([0,1],1,payloadBits_per_OFDM);
array = 3:-1:0;
bin = 2.^(array);
SP = reshape(Txbits,[length(datacarriers),mu]);
index = SP*bin'+1;
%MODULATION
QAM = s(index);
%PLOTTING CONSTELLATION POINTS
figure();
plot(real(s),imag(s),'ro','LineWidth',2)
axis([-4 4 -4 4])
xlabel('Real part')
ylabel('Imaginary part')
title('16 QAM Constellation Diagram')
%Total of 220 bits are being transmitted
%which are converted to 55X4 
%OFDM data with 55 Symbols and 9 Pilots, so Toatal = 64(=K)
ofdm_data = zeros(K,1);
%Loading pilot values at choosen indices
ofdm_data(pilotcarriers) = pilotValue;
%Rest 55 will be data carriers (QAM converterd complex symbols)
ofdm_data(datacarriers)  = QAM;

%CONVERTING OFDM TO TIME DOMAIN
ofdm_time = ifft(ofdm_data);

%CYCLIC PREFIXES
%cyclic prefix will be last 16 values from time-domain ofdm data
cp = ofdm_time(end-15:end)';
%PARALLEL TO SERIAL CONVERTER
ofdm_time =transpose(ofdm_time);
%ADDING CYCLIC PREFIXES;
ofdm_cp = horzcat(cp,ofdm_time);
%TRANSMITTER
%OFDM SYMBOLS ARE TRANSMITTED THROUGH AWGN CHANNEL

%CHANNEL with Impulse response 
channelResponse = [1, 0, 0.3+0.3j, 0.5j]; %the impulse response of the wireless channel
H_exact = fft(channelResponse, K);
figure();
plot(abs(H_exact),'Linewidth',1.5);grid on;
xlabel('Frequency');
ylabel('|H(f)|');
title('Magnitude Spectrum of Channel');
%DATA POWER
snr =30;%signal to noise-ratio in dB at the receiver

Ch_ofdm = conv(ofdm_cp, channelResponse);
signal_pwr = mean(abs(Ch_ofdm.^2));
noise_pwr = signal_pwr.*10.^(-snr/10);
sigma = sqrt(noise_pwr/2);
%AWGN noise
noise = sigma*randn(1,length(Ch_ofdm))+1i*sigma*randn(1,length(Ch_ofdm));
Tx = Ch_ofdm+noise;
figure();
plot(abs(Ch_ofdm),'r-','linewidth',1.5);hold on;
plot(abs(Tx),'b-','linewidth',1.5)
legend('TX signal','RX signal')
title('Transmitted Signal and Recieved Signal at SNR of 20dB')
grid on;
%MEASURE SNR
snr_meas = 10*log10(mean(abs(Ch_ofdm.^2))/mean(abs(noise.^2)));

%RECIEVER DESIGN
%CP is removed 
Rx = Tx(CP+1:(CP+K));
%SERIAL TO PARALLEL CONVERTER
Rx = Rx(:);
%CONVERTED BACK TO FREQUENCY DOMAIN
Rx_fft = fft(Rx);
%ESTIMATING PILOTS AT RECIEVER
pilot = Rx_fft(pilotcarriers);
Rcvd_data = setdiff(Rx_fft,pilot);
Hest_at_pilots = transpose(pilot) / pilotValue;
Hest_abs = interp1(pilotcarriers, abs(Hest_at_pilots), allcarriers);
Hest_phase = interp1(pilotcarriers, angle(Hest_at_pilots), allcarriers);
Hest = Hest_abs .* exp(1i*Hest_phase);  
figure();
plot(allcarriers, abs(H_exact),'linewidth',1.5);
hold on;
stem(pilotcarriers, abs(Hest_at_pilots),'linewidth',1.5);
hold on;
plot(allcarriers, abs(Hest),'m','linewidth',1.5 );
legend('Correct Channel','Pilot estimates','Estimated channel via interpolation');
Eq_hest = Rx_fft ./ transpose(Hest);
Eq_hest(pilotcarriers) =[];
figure();
plot(real(Eq_hest), imag(Eq_hest), 'bo','linewidth',1.5);
hold on;
plot(real(s), imag(s), 'ro','linewidth',2);
axis([-4 4 -4 4])
%Demodulating
e = Eq_hest-s;  
e1 = abs(e).^2;
e1 = transpose(e1);
[p,q] = min(e1);
SER= SER+ sum(q~=index');