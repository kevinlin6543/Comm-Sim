%% ECE 300
clear all;
close all;
clc;
%% Transmission with Viterbi Encoding
%% Setup for Transmission
numIter = 10;  % The number of iterations of the simulation    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

traceBack = 32;
trellis(1) = poly2trellis(7, [171 133]);
trellis(2) = poly2trellis(7, [101 123]);

% The number of symbols per packet of Viterbi encoding
nSym = 1000*trellis(1).numInputSymbols/trellis(1).numOutputSymbols;

M = 4;        % The M-ary number, 2 corresponds to binary modulation
k = log2(M);

% Different Channel that were used in testing the BER
%chan = 1;         
chan = [1 .2 .4]; 
%chan = [0.227 0.460 0.688 0.460 0.227]';  

berVec = zeros(numIter, lenSNR);
berVecQAM = zeros(numIter, lenSNR);
 
%% Simulation
% Generates random bits, reshapes the bits for modulation and then add
% noise and channel. The receiver has an equalizer. Different equalizers
% were shown and the BER was simulated
for i = 1:numIter
   
    bits = randi(2,[nSym*k, 1])-1; 
   
    msg(:,1) = convenc(bits,trellis(1));
    msg(:,2) = convenc(bits,trellis(2));
    for j = 1:lenSNR 
        tx = qammod(msg,M, 'InputType', 'bit','UnitAveragePower',true); 
        
        % Chooses which channel is used
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,tx);
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
       
        % Scale the noise to match for each symbol
        if (M == 2)
            txNoisy = awgn(txChan,3+SNR_Vec(j),'measured'); % Add AWGN
        else 
            txNoisy = awgn(txChan,10*log10(k)+SNR_Vec(j),'measured'); 
        end
        
        rx = qamdemod(txNoisy,M,'OutputType','bit','UnitAveragePower',true);
        
        dataRx(:,1) = vitdec(rx(:,1),trellis(1),traceBack,'cont','hard');
        dataRx(:,2) = vitdec(rx(:,2),trellis(2),traceBack,'cont','hard');
    
        rxMSG = dataRx;
        [~, berVec(i,j,1)] = biterr(bits(1:end-traceBack), rxMSG(traceBack+1:end,1));
        ber(:,1) = mean(berVec(:,:,1));
        [~, berVec(i,j,2)] = biterr(bits(1:end-traceBack), rxMSG(traceBack+1:end,2));
        ber(:,2) = mean(berVec(:,:,2));
        
    end  % End SNR iteration
end      % End numIter iteration
 
%% Plot for BERs
% Takes the mean BER from the demod and constructs graph comparing
% different equalizer and symbol training
berTheory = berawgn(SNR_Vec,'qam',M);
figure();
semilogy(SNR_Vec, ber(:,1:2));
hold on
semilogy(SNR_Vec,berTheory,'r');
title('Viterbi algorithm used on a 16 QAM signal');
xlabel('SNR','fontsize',18);
ylabel('BER','fontsize',18);
legend({'Theoretical 4-QAM', 'Trellis code: 171 133', 'Trellis code: 101 123'});
hold off;

%% Transmission with RS Encoding and 8PSK
clear all;
%% Setup for Transmission
numIter = 10;  % The number of iterations of the simulation
nSym = 700;    % % The number of symbols per packet of RS encoding
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
 
% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

% Create Reed-Solomon Encoder/Decoder
rsEncoder = comm.RSEncoder('BitInput',true,'CodewordLength',N,'MessageLength',K);
rsDecoder = comm.RSDecoder('BitInput',true,'CodewordLength',N,'MessageLength',K);

%% Simulation
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
        else
            txChan = filter(chan,1,tx);  % Apply the channel.
        end
       
        % Convert from EbNo to SNR.
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

%% Plot for BERs
% Takes the mean BER from the demod and constructs graph comparing
% different equalizer and symbol training
berPSK(:,1) = mean(berVec, 1);
figure();
semilogy(SNR_Vec, berPSK);
berTheory = berawgn(SNR_Vec,'psk',M, 'nodiff');
hold on;
semilogy(SNR_Vec,berTheory,'r')
title('Reed Solomon algorithm used on a 8 PSK signal');
xlabel('SNR','fontsize',18);
ylabel('BER','fontsize',18);
legend({'Theoretical BER', 'Actual BER'});
hold off;

%% Transmission with BCH Encoding with 64QAM 
clear all;
%% Setup for Transmission
numIter = 10;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

M = 64;  % The M-ary number
N = 63;  % Codeword length
K = 51;  % Message length
S = 39;  % Shortened message length
k = log2(M);

%chan = 1;          % No channel
chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
%chan = [0.227 0.460 0.688 0.460 0.227]';   % Not so invertible, severe ISI

gp = bchgenpoly(N,K,[]);
bchEncoder = comm.BCHEncoder(N,K,gp,S);
bchDecoder = comm.BCHDecoder(N,K,gp,S);
% Create a vector to store the BER computed during each iteration
berVec = zeros(numIter, lenSNR);

%% Simulation
for i = 1:numIter
    
    %  bits = randint(1, nSym*M, [0 1]);     % Generate random bits
    bits = randi(2,[S*k*nSym, 1])-1;
    
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
        encData = bchEncoder(bits);
        txSig = qammod(encData,M, ...
            'UnitAveragePower',true,'InputType','bit');
        if isequal(chan,1)
            txChan = tx;
        elseif isa(chan,'channel.rayleigh')
            reset(chan) % Draw a different channel each iteration
            txChan = filter(chan,txSig);
        else
            txChan = filter(chan,1,txSig);  % Apply the channel.
        end
        
        % Convert from EbNo to SNR.
        txNoisy = awgn(txChan,10*log10(k)+SNR_Vec(j),'measured'); % Add AWGN
        
        demodSig = qamdemod(txNoisy,M, ...
            'UnitAveragePower',true,'OutputType','bit');
        rxData = bchDecoder(demodSig);
        
        % Compute and store the BER for this iteration
        
        [~, berVec(i,j)] = biterr(bits, rxData);  % We're interested in the BER, which is the 2nd output of BITERR
        
    end  % End SNR iteration
end      % End numIter iteration

% Compute and plot the mean BER
ber = mean(berVec,1);
figure();
semilogy(SNR_Vec, ber)

%% Plot for BERs
% Takes the mean BER from the demod and constructs graph comparing
% different equalizer and symbol training
berTheory = berawgn(SNR_Vec,'qam',M);
hold on
semilogy(SNR_Vec,berTheory,'r')
title('BCH algorithm used on a 64 QAM signal');
xlabel('SNR','fontsize',18);
ylabel('BER','fontsize',18);
legend({'Theoretical BER', 'Actual BER'});
hold off;