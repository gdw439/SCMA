clc;clear;
close all;
% % 同步序列
% FreWidth = 128;
% nzc = exp(1j.*pi./FreWidth.*[0:(FreWidth-1)].^2);                                    % ZC序列
% SynBlock = real(ifft(nzc).').';
% SynPart = [SynBlock(end-16+1:end),SynBlock];                         % 添加CP
% plot(conv(SynPart, SynPart))

% N = 1024;
% chu_ = exp(1j.*pi./N.*[0:(N-1)].^2);
% chu__fft = fft(chu_.').';
% 
% trainseq = repmat(chu_,1,3);
% trainseq_fft = repmat(chu__fft,1,3);
% 
% for i = 1:N
%     correlation(i) = sum(trainseq((N/2+i):(N/2+i+N-1)).*conj(chu_));
% end
% for i = 1:N
%     correlation_fft(i) = sum(trainseq_fft((N/2+i):(N/2+i+N-1)).*conj(chu__fft));
% end
% figure(1)
% stem([1:N],abs(correlation),'r');
% figure(2)
% stem([1:N],abs(correlation_fft),'k')

% 同步序列
LENS = 63;
nzc = exp(1j.*pi./LENS.*(0:(LENS-1)).^2);                                    % ZC序列
SynBlock = ifft([0, nzc, 0, conj(nzc(end:-1:1))]);
SynPart = [SynBlock(end-16+1:end),SynBlock];                         % 添加CP
plot(abs([0, nzc, 0, conj(nzc(end:-1:1))]))