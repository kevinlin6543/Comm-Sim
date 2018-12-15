%% ECE 300
clear all;
close all;
clc;
%% Setup for Transmission
% Determines the number of iterations, number of symbols, the SNR range,
% what modulation is used and the amount of training symbols. The setup
% variables are declared here
numIter = 1000;  % The number of iterations of the simulation
nSym = 1000;    % The number of symbols per packet
SNR_Vec = 0:2:16;
lenSNR = length(SNR_Vec);

traceBack = 32;
trellis = poly2trellis(7,[171 133]);
rate = 1/2;

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
   
    msg = convenc(bits,trellis);
 
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
        
        dataRx = vitdec(rx,trellis,traceBack,'cont','hard');
        
    
        rxMSG = dataRx;
        [~, berVec(i,j,1)] = biterr(bits(1:end-traceBack), rxMSG(traceBack+1:end));
        ber(:,1) = mean(berVec(:,:,1));
        [~, berVec(i,j,2)] = biterr(msg(1:end-traceBack), rx(traceBack+1:end));
        ber(:,2) = mean(berVec(:,:,2));
        
       
    end  % End SNR iteration
end      % End numIter iteration
 
%% Plot for BERs
% Takes the mean BER from the demod and constructs graph comparing
% different equalizer and symbol training
berTheory = berawgn(SNR_Vec,'qam',M);
figure;
semilogy(SNR_Vec, ber(:,1:2));
hold on
semilogy(SNR_Vec,berTheory,'r');
xlabel('SNR')%,'fontsize',18);
ylabel('BER')%,'fontsize',18);