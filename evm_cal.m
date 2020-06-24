% 计算EVM值
% valueR: 接收端信号;
% valueT: 发送端信号;

function evm = evm_cal(valueT,valueR)
    if length(valueR)~=length(valueT)
        error('数据长度不一致！');
    end
    Verr = valueT-valueR;
    Eerr =mean(Verr.*conj(Verr));
    Estd = mean(abs(valueT).^2);
    evm = Eerr/Estd;
    evm = 10*log10(evm);
end