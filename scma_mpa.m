% Saint Petersburg Electrotechnical University, Saint Petersburg, Russia
% Faculty of Radio Engineering
% Department of Theoretical Fundamentals of Radio Engineering
% Vyacheslav P. Klimentyev and Alexander B. Sergienko, 2015

function LLR = scma_mpa2(y, CB, N0, Heff, Nit)
%  SCMA decoder (Log-MPA)
%
%  K - number of orthogonal resources
%  V - number of users (layers)
%  M - number of codewords in each codebook
%
%  N - frame size
%
%  Input arguments:
%
%  y   - SCMA signal after fading channel (size KxN)
%  CB  - SCMA codebooks (size KxMxJ)
%  h   - channel coefficients (size KxJxN)
%  N0  - variance of noise (in AWGN channel)
%  Nit - number of MPA iterations
%
%  Output arguments:
%
%  LLR - Log-Likelihood Ratio (size (log2(M)*V)xN)

map = [0 1 1 0 1 0;
       1 0 2 0 0 1;
       0 2 0 1 0 2;
       2 0 0 2 2 0];
   
K = size(CB, 1);
M = size(CB, 2);
V = size(CB, 3);

% Factor graph calculation
F = zeros(K, V);
s = [K, M];
for k = 1:V
    IND = find(CB(:,:,k));
    [I, ~] = ind2sub(s, IND);
    F(unique(I),k) = 1;
end

debug = [];
N   = size(y, 2);
% LLR = zeros(log2(M)*V, N);

[d_r, d_c] = size(y);
d_c = d_c / 4;

% Step 1: Initial calculations
for jj = 0:1:d_r-1
    for kk = 0:1:d_c-1
        f = zeros(M, M, M, K);
        tmp = y(jj+1, kk*4+1: kk*4+4);
        for k = 1:K % resourses
            ind = find(F(k,:)==1); % non-zero elements, paths
            for m1 = 1:M
                for m2 = 1:M
                    for m3 = 1:M
    %                     f(m1,m2,m3,k) = -(1/N0)*abs(y(k,jj)-(CB(k,m1,ind(1))*h(k,ind(1),jj)+CB(k,m2,ind(2))*h(k,ind(2),jj)+CB(k,m3,ind(3))*h(k,ind(3),jj)))^2;
                        f(m1,m2,m3,k) = -(1/N0)*abs(tmp(k) - (CB(k,m1,ind(1))*Heff(kk*2+map(k, ind(1))+(ind(1)-1)*30) + CB(k,m2,ind(2))* ...
                                        Heff(kk*2+map(k, ind(2))+(ind(2)-1)*30) + CB(k,m3,ind(3))*Heff(kk*2+map(k, ind(3))+(ind(3)-1)*30) ))^2;
%                         f(m1,m2,m3,k) = -(1/N0)*abs(tmp(k)-(CB(k,m1,ind(1))+CB(k,m2,ind(2))+CB(k,m3,ind(3))))^2;
                    end
                end
            end
        end

        Ap = 1/M;
        Igv = zeros(K, V, M); % 资源节点 用户节点 星座图维度 
        Ivg = log(Ap*ones(K, V, M));

    % Step 2: Iterative procedure
        for iter = 1:Nit
            % Igv update
            % 功能节点更新
            for k = 1:K
                % 查找活跃用户
                ind = find(F(k,:)==1);
                    
                for m1 = 1:M
                    sIgv = zeros(1, M*M);
                    for m2 = 1:M
                        for m3 = 1:M
                            sIgv((m2-1)*M+m3) = f(m1,m2,m3,k)+Ivg(k,ind(2),m2)+Ivg(k,ind(3),m3);
                        end
                    end
                    Igv(k,ind(1),m1) = log_sum_exp(sIgv);
                end

                for m2 = 1:M
                    sIgv = zeros(1, M*M);
                    for m1 = 1:M
                        for m3 = 1:M
                            sIgv((m1-1)*M+m3) = f(m1,m2,m3,k)+Ivg(k,ind(1),m1)+Ivg(k,ind(3),m3);
                        end
                    end
                    Igv(k,ind(2),m2) = log_sum_exp(sIgv);
                end

                for m3 = 1:M
                    sIgv = zeros(1, M*M);
                    for m1 = 1:M
                        for m2 = 1:M
                            sIgv((m1-1)*M+m2) = f(m1,m2,m3,k)+Ivg(k,ind(1),m1)+Ivg(k,ind(2),m2);
                        end
                    end
                    Igv(k,ind(3),m3) = log_sum_exp(sIgv);
                end
            end

            % Ivg update 
            % 用户节点更新
            for k = 1:V
                ind = find(F(:,k)==1);
                for n = 1:M
                    Ivg(ind(1),k,n) = log(Ap)+Igv(ind(2),k,n)-log(sum(exp(Igv(ind(2),k,:))));
                    Ivg(ind(2),k,n) = log(Ap)+Igv(ind(1),k,n)-log(sum(exp(Igv(ind(1),k,:))));
                end
            end
        end

    % Step 3: LLR calculation
        Q = zeros(M, V);
        for k = 1:V
            ind = find(F(:,k)==1);
            for m = 1:M
                Q(m,k) = log(Ap)+Igv(ind(1),k,m)+Igv(ind(2),k,m);
            end
        end

        for k = 1:V
            LLR(2*k-1,jj*d_c+kk+1) = log((exp(Q(1,k))+exp(Q(2,k)))/((exp(Q(3,k))+exp(Q(4,k)))));
            LLR(2*k,  jj*d_c+kk+1) = log((exp(Q(1,k))+exp(Q(3,k)))/((exp(Q(2,k))+exp(Q(4,k)))));
        end
    end
end
