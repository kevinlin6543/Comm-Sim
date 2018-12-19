%% ECE 300
clear all;
close all;
clc;
% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
numIter = 10;  % The number of iterations of the simulation
nSym = 700;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);
 
M = 8;  % The M-ary number, we are using 8-PSK
bps = log2(M);
N = 7;  % Codeword length 7
K = 5;  % Message length 5
        % Above numbers result to a total 980 symbols before sending,
        % this is because we can only send 1000 symbols maximum. 980 
        % comes as 700 / 5 * 7


%chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI

 
% Time-varying Rayleigh multipath channel, try it if you dare. Or take
% wireless comms.
% ts = 1/1000;
% chan = rayleighchan(ts,1);
% chan.pathDelays = [0 ts 2*ts];
% chan.AvgPathGaindB = [0 5 10];
% chan.StoreHistory = 1; % Uncomment if you want to be able to do plot(chan)
%
 
% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

% Create Reed-Solomon Encoder/Decoder
rsEncoder = comm.RSEncoder('BitInput',true,'CodewordLength',N,'MessageLength',K);
rsDecoder = comm.RSDecoder('BitInput',true,'CodewordLength',N,'MessageLength',K);

% Run the simulation numIter amount of times
for i = 1:numIter

    bits = randi(2, [nSym*bps, 1])-1;       % Generate random bits
    % New bits must be generated at every
    % iteration
   
    % If you increase the M-ary number, as you most likely will, you'll need to
    % convert the bits to integers. See the BIN2DE function
    % For binary, our MSG signal is simply the bits
    encData = rsEncoder(bits);
    msg = bi2de(reshape(encData, [length(encData)/bps bps]));
    
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
        
        tx = pskmod(msg,M);     % 8-PSK modulate the signal
        
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
       
        % Convert from EbNo to SNR.
        % Note: Because No = 2*noiseVariance^2, we must add ~3 dB
        % to get SNR (because 10*log10(2) ~= 3).
        if (M == 2)
            txNoisy = awgn(txChan,3+SNR_Vec(j),'measured'); % Add AWGN
        else 
            txNoisy = awgn(txChan,10*log10(bps)+SNR_Vec(j),'measured'); % Add AWGN
        end
        
        rx = pskdemod(txNoisy, M);
        rx = de2bi(rx);
        rxMSG = reshape(rx, [length(encData), 1]);
        rxMSG = rsDecoder(rxMSG);
        [~, berVec(i,j)] = biterr(bits, rxMSG);
        
    end  % End SNR iteration
end      % End numIter iteration
 

% Compute and plot the mean BER
berPSK(:,1) = mean(berVec, 1);

figure();
semilogy(SNR_Vec, berPSK);

berTheory = berawgn(SNR_Vec,'psk',M, 'nodiff');
hold on;
semilogy(SNR_Vec,berTheory,'r')
hold off;