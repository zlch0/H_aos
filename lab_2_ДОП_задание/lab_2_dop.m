%% ЭТОТ блок - ЛАБА 2, ниже блоки доп заданий
clc;
clear;
close all;

rng(1);                        

Nfft = 2048;
T = 100e-6;
Fs = Nfft / T;
NumSC = 1600;
NumGuardSC = round((Nfft - NumSC)/2);
all_subc = (NumGuardSC+1):(NumGuardSC+NumSC);
M_data  = 16;                      
M_pilot = 4;                       
SNR_dB = 0:2:20;

NumOfdmSymbols  = 30; % символов на 1 реализацию
numRealizations = 20; % независимых прогонов канала на точку

pilot_step = 8;

kDataBits  = log2(M_data);
kPilotBits = log2(M_pilot);

CP = Nfft/4;

fd_values = [0, 50, 2000]; % Гц

profiles = {
    struct('name','EPA', 'delays',[0 30 70 90 110 190 410]*1e-9, ...
           'gains',[0 -1 -2 -3 -8 -17.2 -20.8]);
    struct('name','EVA', 'delays',[0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
           'gains',[0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]);
    struct('name','ETU', 'delays',[0 50 120 200 230 500 1600 2300 5000]*1e-9, ...
           'gains',[-1 -1 -1 0 0 0 -3 -5 -7]);
};

pilot_idx    = 1:pilot_step:NumSC;
data_idx     = setdiff(1:NumSC, pilot_idx);
nPilots      = length(pilot_idx);
nDataSymbols = length(data_idx);

% суммарное число бит на SNRточку
totalBits = numRealizations * NumOfdmSymbols * nDataSymbols * kDataBits;

BER_res = zeros(length(fd_values), length(SNR_dB), length(profiles));

for idxP = 1:length(profiles)
    prof = profiles{idxP};
    fprintf('\nКанал %s (CP = %d)\n', prof.name, CP);

    for idxF = 1:length(fd_values)
        fd = fd_values(idxF);
        fprintf('  f_d = %d Гц (%d прогонов × %d символов)\n', ...
                fd, numRealizations, NumOfdmSymbols);

        % накопитель ошибок errBits(r, k)
        errBits = zeros(numRealizations, length(SNR_dB));

        % внешний цикл -независимые прогоны 
        for r = 1:numRealizations

            % SNr цикл - каждый (r, k) свой канал и биты 
            for k = 1:length(SNR_dB)

                % новый канал на каждый (r, k), те независимые прогоны
                rayleighchan = comm.RayleighChannel(...
                    'SampleRate',          Fs, ...
                    'PathDelays',          prof.delays, ...
                    'AveragePathGains',    prof.gains, ...
                    'MaximumDopplerShift', fd, ...
                    'NormalizePathGains',  true);

                errSym = 0;

                % один канал на 30 OFDM символов
                for sym = 1:NumOfdmSymbols

                    % данные
                    bits_tx  = randi([0 1], nDataSymbols * kDataBits, 1);
                    mod_data = qammod(bits_tx, M_data, 'UnitAveragePower', true, ...
                                      'InputType', 'bit');

                    % пилоты
                    bits_pilot  = randi([0 1], nPilots * kPilotBits, 1);
                    mod_pilots  = qammod(bits_pilot, M_pilot, 'UnitAveragePower', true, ...
                                         'InputType', 'bit');

                    % OFDMьсимвол 
                    Sp = zeros(Nfft, 1);
                    Sp(all_subc(data_idx))  = mod_data(:);
                    Sp(all_subc(pilot_idx)) = mod_pilots(:);

                    Sig = ifft(ifftshift(Sp));

                    % CP
                    Tx_Signal = [Sig(end-CP+1:end); Sig];

                    %  канал(тотже объект для всех 30 символов) 
                    Rx_Signal = rayleighchan(Tx_Signal);

                    % AWGN
                    Rx_Signal = awgn(Rx_Signal, SNR_dB(k), 'measured');

                    % снятие CP
                    Sig_rec = Rx_Signal(CP+1 : CP+Nfft);

                    % FFT
                    Sp_rec = fftshift(fft(Sig_rec, Nfft));

                    % оценка по пилотам + интерполяция
                    H_pilots = Sp_rec(all_subc(pilot_idx)) ./ mod_pilots(:);
                    H_est = interp1(pilot_idx(:), H_pilots(:), (1:NumSC).', ...
                                    'linear', 'extrap');

                    % ZF
                    Sp_eq = Sp_rec(all_subc(data_idx)) ./ H_est(data_idx);

                    %демодуляция
                    rx_bits = qamdemod(Sp_eq, M_data, 'UnitAveragePower', true, ...
                                       'OutputType', 'bit');

                    errSym = errSym + sum(bits_tx ~= rx_bits);
                end

                errBits(r, k) = errSym;
            end
        end

        % ошибок/бит 
        for k = 1:length(SNR_dB)
            BER_res(idxF, k, idxP) = sum(errBits(:, k)) / totalBits;
        end
    end
end


%% ДОП: п. 1 - ИХ и ЧХ каналов
for idxP = 1:length(profiles)
    prof = profiles{idxP};

    % ИХ
    h_imp = zeros(1, round(max(prof.delays)*Fs) + 2);
    for l = 1:length(prof.delays)
        n = round(prof.delays(l) * Fs) + 1;
        h_imp(n) = h_imp(n) + 10^(prof.gains(l)/20); % лин амплитуда
    end

    % только реальные лучи (без нулевых отсчетов)
    idx = find(h_imp > 0);
    t_imp = (idx - 1) / Fs * 1e6;  %мкс
    amp_dB = 20*log10(abs(h_imp(idx)));%один раз в дб

    figure(1); 
   subplot(3,2,(idxP-1)*2 + 1);
stem(t_imp, h_imp(idx), 'filled', 'LineWidth', 1.5);
grid on; xlabel('задержка, нс'); ylabel('амплитуда');
title([prof.name ': ИХ (лин. масштаб)']); ylim([0 2]);

subplot(3,2,(idxP-1)*2 + 2);
stem(t_imp, 20*log10(h_imp(idx)), 'filled', 'LineWidth', 1.5);
grid on; xlabel('ззадержка, нс'); ylabel('|h|, дБ');
title([prof.name ': ИХ (в дБ)']); ylim([-40 5]);

    % ЧХ
    H = fftshift(fft(h_imp, Nfft));
    H = H / max(abs(H));                       
    f_axis = (-Nfft/2:Nfft/2-1) * (Fs/Nfft) / 1e6;% мгц

    figure(2); subplot(3,1,idxP);
    plot(f_axis, 20*log10(abs(H)+eps), 'LineWidth', 1.2);
    grid on;
    xlabel('f, МГц');
    ylabel('|H|, дБ');
    title(sprintf('%s: АЧХ (fd = 0 Гц)', prof.name));
    xlim([-Fs/2e6 Fs/2e6]);
    ylim([-40 5]);

end

%% ДОП: п. 1 - созвездия
idxF_demo = 1;                 
k_demo    = length(SNR_dB);    
fd_demo   = fd_values(idxF_demo);

for idxP_demo = 1:length(profiles)   
    prof_demo = profiles{idxP_demo};

    %свой объект на каждый профиль
    rayleighchan_demo = comm.RayleighChannel(...
        'SampleRate',          Fs, ...
        'PathDelays',          prof_demo.delays, ...
        'AveragePathGains',    prof_demo.gains, ...
        'MaximumDopplerShift', fd_demo, ...
        'NormalizePathGains',  true);

    % OFDM символ
    bits_tx  = randi([0 1], nDataSymbols * kDataBits, 1);
    mod_data = qammod(bits_tx, M_data, 'UnitAveragePower', true, ...
                      'InputType', 'bit');

    bits_pilot = randi([0 1], nPilots * kPilotBits, 1);
    mod_pilots = qammod(bits_pilot, M_pilot, 'UnitAveragePower', true, ...
                        'InputType', 'bit');

    Sp = zeros(Nfft, 1);
    Sp(all_subc(data_idx))  = mod_data(:);
    Sp(all_subc(pilot_idx)) = mod_pilots(:);

    Sig = ifft(ifftshift(Sp));
    Tx_Signal = [Sig(end-CP+1:end); Sig];

    Rx_Signal = rayleighchan_demo(Tx_Signal);
    Rx_Signal = awgn(Rx_Signal, SNR_dB(k_demo), 'measured');

    Sig_rec = Rx_Signal(CP+1 : CP+Nfft);
    Sp_rec  = fftshift(fft(Sig_rec, Nfft));

    H_pilots = Sp_rec(all_subc(pilot_idx)) ./ mod_pilots(:);
    H_est = interp1(pilot_idx(:), H_pilots(:), (1:NumSC).', ...
                    'linear', 'extrap');

    rx_before = Sp_rec(all_subc(data_idx));
    rx_after  = rx_before ./ H_est(data_idx);

    figure('Name', sprintf('созвездие %s', prof_demo.name), 'Color', 'w');

    subplot(1,2,1);
    plot(real(rx_before), imag(rx_before), '.', 'MarkerSize', 4);
    hold on;
    plot(real(mod_data), imag(mod_data), 'r+', 'MarkerSize', 8, 'LineWidth', 1);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title('до ZF');

    subplot(1,2,2);
    plot(real(rx_after), imag(rx_after), '.', 'MarkerSize', 4);
    hold on;
    plot(real(mod_data), imag(mod_data), 'r+', 'MarkerSize', 8, 'LineWidth', 1);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title('после ZF');

    sgtitle(sprintf('КАМ-%d: %s, f_d = %d Гц, SNR = %d дБ', ...
            M_data, prof_demo.name, fd_demo, SNR_dB(k_demo)));
end
%% ДОП: п. 2 - ber для EPA/EVA/ETU fd=0
idxF_zero = find(fd_values == 0, 1); 

figure();
hold on; grid on;

chan_markers = {'k-o', 'b-s', 'r-d'};
chan_colors  = {'k', 'b', 'r'};

for idxP = 1:length(profiles)
    prof = profiles{idxP};
    Y = reshape(BER_res(idxF_zero, :, idxP), [], 1);

    semilogy(SNR_dB, Y, chan_markers{idxP}, ...
             'Color', chan_colors{idxP}, 'LineWidth', 1.5, 'MarkerSize', 8, ...
             'DisplayName', sprintf('%s (f_d = 0 Гц)', prof.name));
end

xlabel('SNR (дБ)');
ylabel('BER');
title(sprintf('BER EPA/EVA/ETU при f_d = 0 Гц (КАМ-%d, CP = %d)', ...
      M_data, CP));
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-4 1]);
hold off;

%% ДОП: п.3 - доплер спектры EPA5 EVA70 ETU300
fd_spec = [5, 70, 300]; 
names_spec = {'EPA', 'EVA', 'ETU'};

%для PathGains
Fs_gains = Fs / (Nfft + CP); 

figure();

for i = 1:length(names_spec)
    prof = profiles{i};
    fd = fd_spec(i);
    
    chan_dop = comm.RayleighChannel(...
        'SampleRate',          Fs, ...
        'PathDelays',          prof.delays, ...
        'AveragePathGains',    prof.gains, ...
        'MaximumDopplerShift', fd, ...
        'NormalizePathGains',  true, ...
        'PathGainsOutputPort', true);
    
    Nsym = 20000; 
    pathGains_dop = zeros(Nsym, length(prof.delays));
    
    for s = 1:Nsym
        tx_probe = complex(randn(Nfft+CP, 1), randn(Nfft+CP, 1));
        [~, pg] = chan_dop(tx_probe);
        pathGains_dop(s, :) = pg(end, :);
    end
    
    ray1 = pathGains_dop(:, 1); 
    
    subplot(3, 1, i);
    
    [Pxx, F] = pwelch(ray1, [], [], [], Fs_gains, 'centered');
    
    Pxx_dB = 10*log10(Pxx / max(Pxx));
    plot(F, Pxx_dB, 'b-', 'LineWidth', 1.2);
    
    xlim([-1.5*fd, 1.5*fd]);
    ylim([-40, 5]);
    
    title(sprintf('%s, f_d = %d Гц (Fs_{gains} = %.0f Гц)', prof.name, fd, Fs_gains));
    xlabel('f, Гц');
    ylabel('норм. мощность, дБ');
    grid on;
end
%% ДОП: п.4 - BER и Tc
markers = {'k-o', 'b-s', 'r-d', 'g-^', 'm-x'};
colors  = {'k', 'b', 'r', 'g', 'm'};

idxP = 1; % канал
prof = profiles{idxP};

% Tc 
fd_Tc = [50, 2000]; % fd
Tsym  = (Nfft + CP) / Fs;% длительность OFDM-символа с CP в сек

fprintf('\nTc для канала %s:\n', prof.name);
for fd = fd_Tc
    Tc     = 1 / fd;   % с
    Tc_ms  = Tc * 1e3;  % мс
    Tc_sym = Tc / Tsym;% в OFDM-символах
    fprintf('  f_d = %4d гц: Tc = %.3f мс (%.2f символов OFDM)\n', ...
            fd, Tc_ms, Tc_sym);
end

figure;
hold on; grid on;

for idxF = 1:length(fd_values)
    fd = fd_values(idxF);
    Y  = reshape(BER_res(idxF, :, idxP), [], 1);

    semilogy(SNR_dB, Y, markers{idxF}, ...
             'Color', colors{idxF}, 'LineWidth', 1.5, 'MarkerSize', 8, ...
             'DisplayName', sprintf('f_d = %d Гц', fd));
end

xlabel('SNR (дб)');
ylabel('BER');
title(sprintf('OFDM %s (CP = %d, КАМ-%d, %d прогонов на %d симв.)', ...
              prof.name, CP, M_data, numRealizations, NumOfdmSymbols));
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-4 1]);
hold off;

%% ДОП: п.5 - MSE оценки канала от ОСШ (канал EPA)
idxP_mse = 1;
prof_mse = profiles{idxP_mse};
fd_mse   = 50;                          
Tsym     = (Nfft + CP) / Fs;

Nreal_mse = 200;
Nsym_mse  = 20;

MSE_res = zeros(1, length(SNR_dB));

fprintf('\nMSE оценки канала: %s, f_d = %d Гц\n', prof_mse.name, fd_mse);

for k = 1:length(SNR_dB)
    snr = SNR_dB(k);
    mse_acc = 0;
    cnt     = 0;

    for r = 1:Nreal_mse
        ch = comm.RayleighChannel(...
            'SampleRate',          Fs, ...
            'PathDelays',          prof_mse.delays, ...
            'AveragePathGains',    prof_mse.gains, ...
            'MaximumDopplerShift', fd_mse, ...
            'NormalizePathGains',  true);

        for sym = 1:Nsym_mse
            % передача
            bits_tx  = randi([0 1], nDataSymbols * kDataBits, 1);
            mod_data = qammod(bits_tx, M_data, 'UnitAveragePower', true, ...
                              'InputType', 'bit');

            bits_pilot  = randi([0 1], nPilots * kPilotBits, 1);
            mod_pilots  = qammod(bits_pilot, M_pilot, 'UnitAveragePower', true, ...
                                 'InputType', 'bit');

            Sp = zeros(Nfft, 1);
            Sp(all_subc(data_idx))  = mod_data(:);
            Sp(all_subc(pilot_idx)) = mod_pilots(:);

            Sig = ifft(ifftshift(Sp));
            Tx_Signal = [Sig(end-CP+1:end); Sig];

            % канал (без шума)
            Rx_clean = ch(Tx_Signal);

            Sig_clean = Rx_clean(CP+1 : CP+Nfft);
            Sp_clean  = fftshift(fft(Sig_clean, Nfft));
            H_true    = Sp_clean(all_subc) ./ Sp(all_subc);% NumSC*1

            % AWGN
            Rx_Signal = awgn(Rx_clean, snr, 'measured');

            %приём  CP + FFT
            Sig_rec = Rx_Signal(CP+1 : CP+Nfft);
            Sp_rec  = fftshift(fft(Sig_rec, Nfft));

            % оценка по пилотам + лин интерполяция
            H_pilots = Sp_rec(all_subc(pilot_idx)) ./ mod_pilots(:);
            H_est    = interp1(pilot_idx(:), H_pilots(:), (1:NumSC).', ...
                               'linear', 'extrap');

            % MSE на поднесущих данных
            err = H_est(data_idx) - H_true(data_idx);
            mse_acc = mse_acc + mean(abs(err).^2);
            cnt     = cnt + 1;
        end
    end

    MSE_res(k) = mse_acc / cnt;
    fprintf('  SNR = %2d дБ: MSE = %.3e\n', snr, MSE_res(k));
end

%
figure();
semilogy(SNR_dB, MSE_res, 'k-o', 'LineWidth', 1.8, 'MarkerSize', 8);
grid on;
xlabel('SNR (дБ)');
ylabel('MSE оценки канала');
title(sprintf('MSE оценки канала: %s, f_d = %d Гц', ...
      prof_mse.name, fd_mse));