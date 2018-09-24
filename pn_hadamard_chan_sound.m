clc;
clear all;
close all;

FS = 200e6; % Used during simulation of up and down conversion
RATE = FS/2; % Sampling rate used for simulation 
NUMBER_OF_SYMBOLS = 1e6;

UPSAMPLING_FACTOR = FS / RATE; 
DOWNSAMPLING_FACTOR = FS / RATE; 
DECIMATION_FACTOR = FS / RATE;

PN_LEN = 1024;
HADAMARD_LEN = 512;
SPREAD_FACTOR = ceil(PN_LEN/HADAMARD_LEN);
NUM_TX = 4;
NUM_RX = 1;

% Raised Cosine Transmit Filter

FILTER_SPAN = 6;    % Filter span in symbol durations
BETA = 0.1;         % Roll-off factor

rctFilt = comm.RaisedCosineTransmitFilter(...
    'Shape',                  'Normal', ...
    'RolloffFactor',          BETA, ...
    'FilterSpanInSymbols',    FILTER_SPAN, ...
    'OutputSamplesPerSymbol', UPSAMPLING_FACTOR);

% Raised Cosine Receive Filter

FILTER_SPAN = 6;           % Filter span in symbol durations
BETA = 0.1;         % Roll-off factor

rcrFilt = comm.RaisedCosineReceiveFilter(...
    'Shape',                  'Normal', ...
    'RolloffFactor',          BETA, ...
    'FilterSpanInSymbols',    FILTER_SPAN, ...
    'InputSamplesPerSymbol', DOWNSAMPLING_FACTOR, ...
    'DecimationFactor', DECIMATION_FACTOR );


h = comm.MIMOChannel;
h.SampleRate = FS;
h.SpatialCorrelation = false; % Independent channels
h.NumTransmitAntennas = NUM_TX;
h.NumReceiveAntennas = NUM_RX;
% h.TransmitCorrelationMatrix = eye(NUM_TX);
% h.ReceiveCorrelationMatrix = [1 0.3;0.3 1];
h.FadingDistribution = 'Rayleigh';
h.PathDelays = [0,1,2,3,4]*10e-9;
h.NormalizePathGains = true;
h.AveragePathGains = [0,-0.9,-4.9,-8,-20];


pn_seq = step(comm.PNSequence('Polynomial',[10 4 3 1 0],'InitialConditions',[1 0 0 0 0 0 0 0 0 0],'SamplesPerFrame',PN_LEN));
pn_seq1 = step(comm.PNSequence('Polynomial',[10 7 3 1 0],'InitialConditions',[1 0 0 0 0 0 0 0 0 0],'SamplesPerFrame',PN_LEN));
pn_seq2 = step(comm.PNSequence('Polynomial',[10 9 8 7 6 5 4 3 0],'InitialConditions',[1 0 0 0 0 0 0 0 0 0],'SamplesPerFrame',PN_LEN));
pn_seq = (2*[pn_seq] - 1);
pn_seq1 = (2*[pn_seq1] - 1);
pn_seq2 = (2*[pn_seq2] - 1);

H_d = hadamard(HADAMARD_LEN);

sig_to_transmit = zeros((PN_LEN+10)*(NUM_TX+1) + NUMBER_OF_SYMBOLS,NUM_TX);

for u = 1:HADAMARD_LEN
    x(:,:,u) = repmat(H_d(u,:),SPREAD_FACTOR,1);
    x_reshaped(:,:,u) = reshape(x(:,:,u),1,SPREAD_FACTOR*HADAMARD_LEN);
    sig(:,u) = transpose(x_reshaped(:,1:PN_LEN,u)) .* pn_seq;
end
[sig_sorted,I] = sort(abs(sum(sig,1)));

K = 4;
for t = 1:NUM_TX
%     x(:,:,t) = repmat(H_d(t,:),SPREAD_FACTOR,1);
%     x_reshaped(:,:,t) = reshape(x(:,:,t),1,SPREAD_FACTOR*HADAMARD_LEN);
%     sig(:,t) = transpose(x_reshaped(:,1:PN_LEN,t)) .* pn_seq;
    X_Input = randi([0 K-1],NUMBER_OF_SYMBOLS,1);
    QAMOutput(:,t) = qammod(X_Input, K);
    sig_to_transmit(10*(t-1)+PN_LEN*(t-1)+1:10*(t-1)+PN_LEN*t,t) = sig(:,I(t));
    sig_to_transmit(((PN_LEN+10)*NUM_TX)+1:((PN_LEN+10)*NUM_TX)+PN_LEN,t) = sig(:,I(t));
    sig_to_transmit(((PN_LEN+10)*NUM_TX)+PN_LEN+11:end,t) = QAMOutput(:,t);
end

y_rct = rctFilt(sig_to_transmit);
y_filt = step(h,y_rct);
y_rcr = rcrFilt([y_filt]);

snr = 3;
y_rcr_noise = awgn(y_rcr,snr);

for r = 1:NUM_RX
%     iter = 0;
%     while iter <= length(y_rcr_noise(:,r))-PN_LEN
%         temp = transpose(y_rcr_noise(1+iter:iter+PN_LEN,r))*pn_seq;
%         out(iter+1) = temp/(PN_LEN);
%         iter = iter+1;
%     end
    for t = 1:NUM_TX
        [rec_corr(:,r,t),lags(:,r,t)] = xcorr(y_rcr_noise(((PN_LEN+10)*NUM_TX)+1:((PN_LEN+10)*NUM_TX)+PN_LEN,r),sig(:,I(t)));
        chan_est(:,r,t) = rec_corr(find(lags(:,r,t)>=0),r,t)/(sig(:,I(t))'*sig(:,I(t)));
        
        THRESHOLD_FOR_CORR_PEAKS = 0.2;

        for i = 1:length(chan_est)
            if abs(chan_est(i,r,t)) < THRESHOLD_FOR_CORR_PEAKS
                chan_est(i,r,t) = 0;
            end
        end
    end
end
for r = 1:NUM_RX
    figure;
    for t = 1:NUM_TX
        plot(abs(chan_est(:,r,t))); hold on;
    end
    title(['RX ',int2str(r)]);
end
for r = 1:NUM_RX
    figure;
    for t = 1:NUM_TX
        plot(lags(:,r,t),abs(rec_corr(:,r,t))./((sig(:,I(t))'*sig(:,I(t))))); hold on;
    end
    title(['RX ',int2str(r)]);
end

