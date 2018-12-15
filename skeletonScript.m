% A skeleton BER script for a wireless link simulation
clear all;close all;clc
% For the final version of this project, you must use these 3
% parameter. You will likely want to set numIter to 1 while you debug your
% link, and then increase it to get an average BER.
numIter = 10;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

berEstSoft = zeros(size(SNR_Vec)); 
berEstHard = zeros(size(SNR_Vec));

trellis = poly2trellis(7,[171 133]);

tbl = 18;
rate = 1/2;

M = 4;        % The M-ary number, 2 corresponds to binary modulation
k = log2(M);

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

% Run the simulation numIter amount of times
for i = 1:numIter
    
    %  bits = randint(1, nSym*M, [0 1]);     % Generate random bits
    bits = randi(2,[nSym*k, 1])-1; 
    % New bits must be generated at every
    % iteration
    
    % If you increase the M-ary number, as you most likely will, you'll need to
    % convert the bits to integers. See the BIN2DE function
    % For binary, our MSG signal is simply the bits
    msg = convenc(bits,trellis);
    
    for j = 1:lenSNR % one iteration of the simulation at each SNR Value
        [numErrsSoft,numErrsHard,numBits] = deal(0);     
        tx = qammod(msg,M, 'InputType', 'bit','UnitAveragePower',true);  % BPSK modulate the signal
        
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
        
        %rx = qamdemod(txNoisy,M); % Demodulate
        
        rxDataHard = qamdemod(txNoisy,M,'OutputType','bit','UnitAveragePower',true);
       
        dataHard = vitdec(rxDataHard,trellis,tbl,'cont','hard');
        
        % Calculate the number of bit errors in the frame. Adjust for the
        % decoding delay, which is equal to the traceback depth.
        [numErrsInFrameHard, berVec(i,j,1)] = biterr(bits(1:end-tbl),dataHard(tbl+1:end));
        %[numErrsInFrameSoft, berVec(i,j,2)] = biterr(bits(1:end-tbl),dataSoft(tbl+1:end));
        ber(:,1) = mean(berVec(:,:,1));
        %ber(:,2) = mean(berVec(:,:,2));
        numErrsHard = numErrsHard + numErrsInFrameHard;
        %numErrsSoft = numErrsSoft + numErrsInFrameSoft;
        numBits = numBits + nSym*k;
        
        % Again, if M was a larger number, I'd need to convert my symbols
        % back to bits here.
        %rxMSG = rx;
        
        % Compute and store the BER for this iteration
        
        %[zzz berVec(i,j)] = biterr(msg, rxMSG);  % We're interested in the BER, which is the 2nd output of BITERR
        
    end  % End SNR iteration
end      % End numIter iteration


% Compute and plot the mean BER
%ber = mean(berVec,1);
semilogy(SNR_Vec, ber(:,1))
% Compute the theoretical BER for this scenario
% THIS IS ONLY VALID FOR BPSK!
% YOU NEED TO CHANGE THE CALL TO BERAWGN FOR DIFF MOD TYPES
% Also note - there is no theoretical BER when you have a multipath channel
berTheory = berawgn(SNR_Vec,'qam',4);
hold on
semilogy(SNR_Vec,berTheory,'r')
legend('BER', 'Theoretical BER')
hold off