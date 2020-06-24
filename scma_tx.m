% parameters: 
%     carriers sum :  63 in 128 carriers
%     DCO-OFDM
%     without constellation rotation, so it is better than the origin .m
%     file.
% by vortex
% 2019.10.21

clc;clear;
close all;
tic

N = 15000;                                                                 % SCMA signals in frame
Fs = 50e6;                                                                 % 采样频率
Niter = 10;
OveSamT = 4; 
Interval = 16;
carrier = 128;
block_len = 100;
block_sum = N / block_len / 15;  

% Codebooks
CB(:,:,1) = [ 0                  0                  0                  0;...
             -0.1815-1j*0.1318  -0.6351-1j*0.4615   0.6351+1j*0.4615   0.1815+1j*0.1318;...
              0                  0                  0                  0;...
              0.7851            -0.2243             0.2243            -0.7851 ];

CB(:,:,2) = [ 0.7851            -0.2243             0.2243            -0.7851;...
              0                  0                  0                  0;...
             -0.1815-1j*0.1318  -0.6351-1j*0.4615   0.6351+1j*0.4615   0.1815+1j*0.1318;...
              0                  0                  0                  0 ];

CB(:,:,3) = [-0.6351+1j*0.4615   0.1815-1j*0.1318  -0.1815+1j*0.1318   0.6351-1j*0.4615;...
              0.1392-1j*0.1759   0.4873-1j*0.6156  -0.4873+1j*0.6156  -0.1392+1j*0.1759;...
              0                  0                  0                  0;...
              0                  0                  0                  0 ];

CB(:,:,4) = [ 0                  0                  0                  0;...
              0                  0                  0                  0;...
              0.7851            -0.2243             0.2243            -0.7851;...
             -0.0055-1j*0.2242  -0.0193-1j*0.7848   0.0193+1j*0.7848   0.0055+1j*0.2242 ];

CB(:,:,5) = [-0.0055-1j*0.2242  -0.0193-1j*0.7848   0.0193+1j*0.7848   0.0055+1j*0.2242;...
              0                  0                  0                  0;...
              0                  0                  0                  0;...
             -0.6351+1j*0.4615   0.1815-1j*0.1318  -0.1815+1j*0.1318   0.6351-1j*0.4615 ];

CB(:,:,6) = [ 0                  0                  0                  0;...
              0.7851            -0.2243             0.2243            -0.7851;...
              0.1392-1j*0.1759   0.4873-1j*0.6156  -0.4873+1j*0.6156  -0.1392+1j*0.1759;...
              0                  0                  0                  0 ];

% Synchronization sequence
Chu = exp(1j.*pi./carrier.*[0:(carrier-1)].^2);                            % ZC序列
SynBlock = real(ifft(Chu));
SynPart = [SynBlock(end-Interval+1:end),SynBlock];                         % 添加CP
          
% Sparse Pilot
SChu       = exp(1j.*pi./60.*(0:1:(60-1)).^2);                             % ZC序列
PilotMask  = repmat(reshape(CB(:, 1, :),4,6).',1, 15) ~= 0;

% 构建测量矩阵
Seff = zeros(60,30,6);
for i = 1:6
    tmp = diag(PilotMask(i,:));
    Seff(:, :, i) = tmp(:, PilotMask(i,:)~=0);
end

S = [Seff(:, :, 1), Seff(:, :, 2), Seff(:, :, 3), Seff(:, :, 4), Seff(:, :, 5), Seff(:, :, 6)];

[~,k] = min(~S,[],2);
tmp = size(S);
I = sub2ind(tmp, 1:tmp(1), k.');
S1 = zeros(tmp);
S1(I) = 1;
S     = S - S1;
S1    = diag(SChu)*S1;

[~,k] = min(~S,[],2);
I = sub2ind(tmp, 1:tmp(1), k.');
S2 = zeros(tmp);
S2(I) = 1;
S3    = S - S2;
S2    = diag(SChu)*S2;
S3    = diag(SChu)*S3;

T = S1'*S1 + S2'*S2 + S3'*S3;

% 构建导频 分配序列
SPilotBox = zeros(3, 60, 6);
for m = 0:5
    SPilotBox(1, :, m+1) = sum(S1(:, 1+m*30:30+m*30), 2).';
    SPilotBox(2, :, m+1) = sum(S2(:, 1+m*30:30+m*30), 2).';
    SPilotBox(3, :, m+1) = sum(S3(:, 1+m*30:30+m*30), 2).';
end

% product the training sequence
Seff = zeros(60,30,6);
for i = 1:6
    tmp1 = SChu.*PilotMask(i,:);
    tmp2 = diag(tmp1(1:60));
    Seff(:, :, i) = tmp2(:, PilotMask(i,:)~=0);
end
    
K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)

% log2(M)-bit symbols
Mbit = randi([0 M-1], V, N);  
save('./Txdata/Mdata.mat', 'Mbit', 'CB', 'SynPart', 'S1', 'S2', 'S3', 'T');

for ue = 1:6     
    word = reshape(CB(:, Mbit(ue, :)+1, ue), 60, []).';
    
    [w_r, w_c]  = size(word);
    
    word_pilot = zeros(block_sum*3+w_r, w_c);
        
     % 导频及前缀
    for block_index = 0 : block_sum - 1
        step_n = block_index*block_len;
        step_o = block_index*(block_len+3);
        word_pilot(step_o + 1 : step_o + 3, :) = SPilotBox(:, :, ue); % SubPilot(ue,:);
        word_pilot(step_o + 4 : step_o + block_len + 3, :) = word(step_n + 1: step_n + block_len, :);
    end
    
    word_conj = conj(word_pilot);

    ofdm_symbol = [zeros(length(word_pilot), 1), word_pilot, zeros(length(word_pilot), 7), word_conj(:, end:-1:1)]; % Hermitian conjugate

    ofdm_signal = ifft(ofdm_symbol, carrier, 2);

    % 添加保护间隔 
    ofdm_interval = [ofdm_signal(:,carrier-Interval+1:end),ofdm_signal];
    % 添加同步头
    ofdm_transmit = [SynPart, reshape(ofdm_interval.', 1,[])];
    % 去直流分量
    ofdm_transmit = ofdm_transmit - mean(ofdm_transmit);
    % AWG只能读列向量
    ofdm_transmit = ofdm_transmit.';
    
    filename = sprintf("./Txdata/ue%d.txt",ue);
    save(filename, 'ofdm_transmit','-ascii');
end
toc