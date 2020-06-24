function [ss_out, tt, Fs_new_out] = oversamp(ss, Fs, Fs_new, Sign_exact_OS, ...
                                f_space, N_samp_add, Monitor)
        % 
        %oversamp Do oversampling or downsampling for received sequency in frequency domain. 
        % 
        %   oversamp(ss, Fs, Fs_new, Sign_OsNotDs) returns oversampled (downsampled)
        %      signal sequency at the frequency of Fs_new. It is implemented by adding
        %      zeros (removing points outside of signal bandwidth) in the frequency
        %
        %   Input: 
        %      ss: received digital signal sequence
        %      Fs: sampling rate of input signal
        %      Fs_new: sampling rate of output signal
        %      Sign_exact_OS:  1: exactly oversample or downsample
        %                      0: roughly ~~~
        %      f_space: frequency space in frequency domain, should be positive
        %         real value.
        %      N_samp_add: the number of samples to be added/substracted into the
        %         sequence, must be even.
        %      Monitor: number of displays for monitoring.
        %               0 -- no display; -1 -- display all.
        %
        %   Output:
        %      ss_out: output signal sequency.
        %      tt: corresponding time sequency.
        %      Fs_new_out: exact new sampling rate of output oversampled signal.

        %   &Revision: 1.2.1 &   & Data: 2008/07/04 &
        %   M-file function

        Sign_addPt = 1;     % Initialization
        if nargin < 7
            Monitor = 0;
            if nargin < 6
                Sign_addPt = 0;
                if nargin < 4
                    Sign_exact_OS = 1;
                    if nargin < 3
                        error('MATLAB:oversamp:NotEnoughInputs',...
                            'Not enough input arguments. See oversamp.')
                    end
                end
            end
        end

        % No oversampling and downsampling, return.
        flag_tmp = 1;
        if Sign_addPt == 1
            if N_samp_add == 0
                flag_tmp = 0;
            end
        else
            if Fs == Fs_new
                flag_tmp = 0;
            end
        end
        if flag_tmp == 0
            ss_out = ss;
            tt = (0: length(ss) - 1)/Fs;
            Fs_new_out = Fs;
            if Monitor >= 1 || Monitor == -1
                disp('No oversampling/downsampling!')
            end
            return
        end

        if f_space < 0 || imag(f_space) ~= 0    % check frequency space value
            error('MATLAB:oversamp:PositiveValue',...
                'f_space should be positive real value. See oversamp.')
        end

        % Initialization
        N_ori = length(ss);
        flag_save = 0;

        if Sign_addPt == 1
           N_new = N_ori + N_samp_add; 
        else
           if Sign_exact_OS == 1      % exactly oversample or downsample
                N_new = (N_ori/Fs)*Fs_new;
            elseif Sign_exact_OS == 0
                N_1 = 2*round(Fs / f_space /2);
        %         N_1 = 2^(nextpow2(N_1));
                if N_1 > N_ori
                    ss_save = ss;
                    N_ori_save = N_ori;
                    ss = [ss zeros(1, N_1 - N_ori)];
                    N_ori = N_1;
                    flag_save = 1;
                end
                if Fs_new >= Fs % oversampling
                    N_new = 2*ceil((N_ori/Fs)*Fs_new/2);  
                else            % downsampling
                    N_new = 2*floor((N_ori/Fs)*Fs_new/2); 
                end
           else
                error('Incorrect parameter setting for if...elseif... command!')
            end
        end
        if mod(N_new, 2) ~= 0
            error('Decimal fraction. N_new should be even!')
        end

        % Oversampling (downsampling) in the frequency domain
        ss_f = fft(ss);
        if N_new >= N_ori    % means Fs_new >= Fs -- oversampling
            ss_f_2 = [ss_f(1: end/2) zeros(1, N_new - N_ori) ss_f(end/2 + 1: end)];
        else                % downsampling
            ss_f_2 = [ss_f(1: N_new/2) ss_f(end - N_new/2 + 1: end)];
        end
        ss_out = (N_new/N_ori) * ifft(ss_f_2);
        Fs_new_out = N_new * (Fs/N_ori);
        tt = (0: N_new - 1)/Fs_new_out; 

        if flag_save == 1
            ss_out = ss_out(1: (N_ori_save + N_new - N_ori));
            tt = tt(1: (N_ori_save + N_new - N_ori));
        end
        if Monitor >= 1 || Monitor == -1
            disp([sprintf('Input sampling rate (Hz): \t\t\t%0.2e\nOutput sampling rate (Hz): \t\t\t%0.2e',Fs, Fs_new_out)])
            tmp = ((1/Fs_new_out) - (1/Fs))*Fs;
            disp([sprintf('Sampling rate difference (Hz): \t\t%0.2e\nSampling clock offset ratio: \t\t%0.2e\n', Fs_new_out - Fs, tmp)])
        end
    end % function [ss_out, tt, Fs_new_out] = oversamp() ----------------------