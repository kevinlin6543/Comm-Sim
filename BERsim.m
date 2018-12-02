%% ECE 300
clear all;
close all;
clc;
% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
numIter = 1;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);
symbolTrain = 200;

M = 16;        % The M-ary number, 2 corresponds to binary modulation
k = log2(M);
%SNR_Vec = SNR_Vec + 10*log10(k);

chan = 1;          % No channel
% chan = [1 .2 .4]; % Somewhat invertible channel impulse response, Moderate ISI
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

% Run the simulation numIter amount of times
for i = 1:numIter
    
    % bits = randi(1, nSym*M, [0 1]);     % Generate random bits
    bits = randi(2,[nSym*k, 1])-1;
    % New bits must be generated at every
    % iteration
    
    % If you increase the M-ary number, as you most likely will, you'll need to
    % convert the bits to integers. See the BIN2DE function
    % For binary, our MSG signal is simply the bits
    msg = bi2de(reshape(bits, [nSym k]));
    %msg(1:50) = mod(1:50, 16);

    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
                
        tx = qammod(msg,M);  % BPSK modulate the signal
        
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
        txNoisy = awgn(txChan,10*log10(k)+SNR_Vec(j),'measured'); % Add AWGN
  
        eq1 = lineareq(2, rls(1));
        eq1.SigConst = qammod((0:M-1)',M)';
        [symbolest, yd] = equalize(eq1, txNoisy, tx(1:symbolTrain));
        eql.ResetBeforeFiltering = 0;
        
        rx = qamdemod(yd,M); % Demodulate
        rx = de2bi(rx); 
        rxMSG = reshape(rx, [nSym*k, 1]);
        
        % Again, if M was a larger number, I'd need to convert my symbols
        % back to bits here.
        % rxMSG = rx;
        
        % Compute and store the BER for this iteration
        
        [~, berVec(i,j)] = biterr(bits(symbolTrain+1:end), rxMSG(symbolTrain+1:end));  % We're interested in the BER, which is the 2nd output of BITERR
        
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean BER
ber = mean(berVec,1);

semilogy(SNR_Vec, ber)

% Compute the theoretical BER for this scenario
% THIS IS ONLY VALID FOR BPSK!
% YOU NEED TO CHANGE THE CALL TO BERAWGN FOR DIFF MOD TYPES
% Also note - there is no theoretical BER when you have a multipath channel
berTheory = berawgn(SNR_Vec,'qam',16);
hold on
semilogy(SNR_Vec,berTheory,'r')
legend('BER', 'Theoretical BER')