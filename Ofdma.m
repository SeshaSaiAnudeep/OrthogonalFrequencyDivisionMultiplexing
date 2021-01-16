clearvars;clc;close all;
% No.of SubCarriers: 64
% CYCLIC PREFIX 
    %The length of the cyclic prefix (CP) denotes the number of samples 
    %that are copied from the end of the modulated block to the beginning
% Modulation: 16-QAM
% NO.OF PILOTS : 9
    %The number of pilots P in the OFDM symbol describes, how many 
    %carriers are used to transmit known information (i.e. pilots).
clearvars;clc;close all;
%16QAM MODULATION
M=16;
k = log2(M);
c= [-3 -1 1 3];
[a b] = meshgrid(c);
s = [reshape(a,[1,M]);reshape(b,[1,M])];
s = s(1,:)+sqrt(-1)*s(2,:);
K = 64; % number of ofdm sub carriers
N=64;%N point fft = #fdm carriers
CP = floor(K/4);%Length of the cyclic prefix: 25% of the block
P = 8;%number of pilor carriers per OFDM block
pilotValue = 3+3i; %the known value of each pilot transmits
% ------------------------------------------ % 
% ###### LTE PHYSICAL LAYER PARAMETERS ##### % 
% ------------------------------------------ %
ofdmBW = 20*10^6; % OFDM bandwidth
deltaF = ofdmBW/N; %=20MHz/64 = 0.3125 MHz
Tfft = 1/deltaF; % IFFT/FFT period = 3.2us
Tgi = Tfft/4;%Guard interval duration - duration of cyclic prefix = 0.8us
Tsignal = Tgi+Tfft; %duration of OFDM symbol =4us
Ncp = N*Tgi/Tfft;%#of symbols allocated to cyclic prefix =16
Nsd = 55;
Nsp = 9;
Nst = Nsd + Nsp; %Number of total used subcarriers
allcarriers = zeros(1,64);
for i =1:64
    allcarriers(i) = i;
end
%===Adding pilot value at every 8th sample===%    
pilotcarriers = allcarriers(1:floor(K/P):end);
% ===== Let the last sample be pilot ===== % 
pilotcarriers = [pilotcarriers allcarriers(end)];
%Remaining indices are deddicated to datacarriers 
datacarriers = setdiff(allcarriers,pilotcarriers);
%16 QAM therefore 4 bits are required to represent each symbol
mu = 4;%BITS PER SYMBOL
payloadBits_per_OFDM = length(datacarriers)*mu; %number of payload bits per OFDM symbol
EsNodB =-20:8;%signal to noise-ratio in dB at the receiver
SER = zeros(1,length(EsNodB));
Nsim = 1000;
Txbits = randi([0,1],1,payloadBits_per_OFDM);
for i = 1:length(EsNodB)
    for sim = 1:Nsim
array = 3:-1:0;
bin = 2.^(array);
SP = reshape(Txbits,[length(datacarriers),mu]);
index = SP*bin'+1;
%MODULATION
QAM = s(index);
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
channelResponse = [1, 0, 0.3+0.3j]; %the impulse response of the wireless channel
H_exact = fft(channelResponse, K);
Ch_ofdm = conv(ofdm_cp, channelResponse);
signal_pwr = mean(abs(Ch_ofdm.^2));
noise_pwr = signal_pwr.*10.^(-EsNodB/10);
sigma = sqrt(noise_pwr/2);
noise = sigma(i)*randn(1,length(Ch_ofdm))+1i*sigma(i)*randn(1,length(Ch_ofdm));
Tx = Ch_ofdm + noise;
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
Eq_hest = Rx_fft ./ transpose(Hest);
Eq_hest(pilotcarriers) =[];
%Demodulating
e = Eq_hest-s;  
e1 = abs(e).^2;
e1 = transpose(e1);
[p,q] = min(e1);
errors = sum(q~=index');
SER(i)= SER(i)+errors;
    end
end
SER = SER/(55*Nsim);
% SER1 = SER/Nsim;
EsNoLin = 10.^(EsNodB/10);
SERth=2*(1-1/sqrt(M))*qfunc(sqrt(3*EsNoLin/(M-1)));
SERth=1-(1-SERth).^2; 
semilogy(EsNodB,SER,'bo-','linewidth',2);
hold on;
semilogy(EsNodB,SERth,'r-','linewidth',2);
grid on;
title('BER')
xlabel('SNRdB')
ylabel('SER')
legend('simulated','theoretical')