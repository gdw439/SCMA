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

Niter = 10;
load ./Txdata/Mdata.mat

K = size(CB, 1); % number of orthogonal resources
M = size(CB, 2); % number of codewords in each codebook
V = size(CB, 3); % number of users (layers)

F = zeros(K, V);
s = [K, M];
for k = 1:V
    IND = find(CB(:,:,k));
    [I, ~] = ind2sub(s, IND);
    F(unique(I),k) = 1;
end

% load user data
ue1 = importdata('./Txdata/ue1.txt').';
ue2 = importdata('./Txdata/ue2.txt').';
ue3 = importdata('./Txdata/ue3.txt').';
ue4 = importdata('./Txdata/ue4.txt').';
ue5 = importdata('./Txdata/ue5.txt').';
ue6 = importdata('./Txdata/ue6.txt').';


N = 15000;                                                                 % SCMA signals in frame
Fs = 50e6;                                                                 % 采样频率
Niter = 6;
OveSamT = 4; 
Interval = 16;
carrier = 128;
block_len = 100;
block_sum = N / block_len / 15;  

rec_signal = ue1 + 0.5*ue2 + ue3 + ue4 + ue5;% +ue6;

syn_signal = synchronize(rec_signal, SynPart, length(SynPart), Interval, block_sum*(block_len+3));

FrameLine   = syn_signal;
FrameRec_GI = reshape(FrameLine, carrier+Interval, []);
FrameRecGT  = FrameRec_GI.';
FrameRec    = FrameRecGT(:,Interval+1:end);                                % 去CP
RPilot_mat1 = FrameRec(1:block_len+3:end,:);                               % 提导频
FrameRec(1:block_len+3:end,:) = [];                                        % 去导频

RPilot_mat2 = FrameRec(1:block_len+2:end,:);                               % 提导频
FrameRec(1:block_len+2:end,:) = [];                                        % 去导频

RPilot_mat3 = FrameRec(1:block_len+1:end,:);                               % 提导频
FrameRec(1:block_len+1:end,:) = [];                                        % 去导频

RPilot_ave = [mean(RPilot_mat1); mean(RPilot_mat2);mean(RPilot_mat3)];
RPilot     = fft(RPilot_ave, carrier, 2);
RSPilot    = RPilot(:, 2:61);

FrameRec = fft(FrameRec, carrier, 2);

Heff = inv(T) * (S1'*RSPilot(1,:).'+ S2'*RSPilot(2,:).' + S3'*RSPilot(3,:).');

TimeSlotCor = FrameRec;
ExtractMap = TimeSlotCor(:, 2 : 1+60);
word = reshape(ExtractMap', 4, []);

y = word;

N0 = 0.001;                                                                % 对于这个参数的估计比较重要
LLR = scma_mpa(ExtractMap, CB, N0, Heff, Niter);% Heff, 
% LLR = scmadec(y, CB, N0, Niter);
% symbol to bit conversion
r    = de2bi(Mbit, log2(M), 'left-msb');
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

err = sum(xor(data, datar));
Nerr  = err.';
Nbits = log2(M)*N;

BER = Nerr./Nbits;

toc