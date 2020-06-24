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


          
K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)

N = 10000; % SCMA signals in frame

EbN0 = 1;
SNR  = EbN0 + 10*log10(log2(M)*V/K);

Nerr  = zeros(V, length(SNR));
Nbits = zeros(V, length(SNR));
BER   = zeros(V, length(SNR));

maxNumErrs = 100;
maxNumBits = 1e6;
Niter      = 10;

carrier = 128;
block_len = 50;
block_sum = N / block_len;  
Fs = 50e6;                                                                 % 采样频率
OveSamT = 4; 
Interval = 16;

FreWidth = 128;
ZC = exp(1j.*pi./FreWidth.*[0:(FreWidth-1)].^2);                           % ZC序列
SynBlock = real(ifft(ZC));
SynPart = [SynBlock(end-Interval+1:end),SynBlock];                         % 添加CP

% 导频
SubCarrier = 3;
n = 1:1:SubCarrier;
nzc = exp( -1i*pi*25*(n-1).*n./SubCarrier);
Nzc = [nzc,nzc(1)];
PilotSeq = zeros(1,carrier);
PilotSeq(2:1+length(Nzc)) = Nzc;
NzcJ = conj(Nzc);
PilotSeq(end:-1:end-length(Nzc)+1) = NzcJ(1:end);   % 共轭
PilotBlock = ifft(PilotSeq);

for k = 1:length(SNR)
    
    while ((min(Nerr(:,k)) < maxNumErrs) && (Nbits(1,k) < maxNumBits))
        N0 = 1/(10^(SNR(k)/10)); % noise power

        x = randi([0 M-1], V, N); % log2(M)-bit symbols

        s = scmaenc(x, CB)'; % joint encoding and fading channel propagation

        s_conj = conj(s);

        ofdm_symbol = [zeros(N, 1), s, zeros(N, carrier-1-K-K), s_conj(:, end:-1:1)]; % Hermitian conjugate

        ofdm_signal = ifft(ofdm_symbol,carrier,2);

        ofdm_pilots = zeros(N + block_sum, carrier); 

         % 导频及前缀
        for block_index = 0 : block_sum - 1
            ofdm_pilots(block_index*(block_len+1) + 1, :) = PilotBlock;
            ofdm_pilots(block_index*(block_len+1) + 2: (block_index+1)*(block_len+1), :) = ofdm_signal(block_index*(block_len) + 1: (block_index+1)*(block_len), :);
        end

        % 保护间隔 
        ofdm_interval = [ofdm_pilots(:,carrier-Interval+1:end),ofdm_pilots];

        ofdm_transmit = [SynPart, reshape(ofdm_interval', 1,[])];

        ofdm_transmit = ofdm_transmit - mean(ofdm_transmit);

        ofdm_recv = awgn(ofdm_transmit, SNR(k), 'measured');

        ofdm_recv = ofdm_recv - mean(ofdm_recv);

        FrameLine = ofdm_recv(145:1:end);
        FrameRec_GI = reshape(FrameLine, carrier+Interval, N+block_sum);
        FrameRecGT = FrameRec_GI.';
        FrameRec = FrameRecGT(:,Interval+1:end);                           % 去CP
        Pilot = FrameRec(1:block_len+1:end,:);                             % 提导频

        FrameRec(1:block_len+1:end,:) = [];                                % 去导频
        FrameRec = fft(FrameRec,carrier,2);
        OriPilot= fft(PilotBlock);
        PilotAve = mean(Pilot);

        H_LS = fft(PilotAve)./OriPilot;
        TimeSlotCor = zeros(N, length(H_LS));
        for floor = 1 : N
            TimeSlotCor(floor,:) = FrameRec(floor,:)./H_LS;
        end
        TimeSlotCor = FrameRec;
        ExtractMap = TimeSlotCor(:, 2 : 1+K);  

        y = ExtractMap';

        LLR = scmadec(y, CB, N0, Niter);

        % symbol to bit conversion
        r    = de2bi(x, log2(M), 'left-msb');
        data = zeros(log2(M)*N, V);
        for kk = 1:V
            data(:,kk) = reshape(downsample(r, V, kk-1).',[],1);
        end

        % LLR to bit conversion
        datadec = reshape((LLR <= 0), [log2(M) N*V]).';
        datar   = zeros(log2(M)*N, V);
        for kk = 1:V
            datar(:,kk) = reshape(downsample(datadec, V, kk-1).', [], 1);
        end

        err        = sum(xor(data, datar));
        Nerr(:,k)  = Nerr(:,k) + err.';
        Nbits(:,k) = Nbits(:,k) + log2(M)*N;
    end
    
    BER(:,k) = Nerr(:,k)./Nbits(:,k);
    k
end

% plot result curve
figure()
semilogy(EbN0,BER(1,:))

toc