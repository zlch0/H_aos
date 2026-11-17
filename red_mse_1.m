clc;
clear;
close all;
%%
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
numRealizations = 100; % независимых прогонов канала на точку

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

% суммарное число бит на SNR-точку
totalBits = numRealizations * NumOfdmSymbols * nDataSymbols * kDataBits;

BER_res = zeros(length(fd_values), length(SNR_dB), length(profiles));
MSE_res = zeros(length(fd_values), length(SNR_dB), length(profiles));

Eb = mean(abs(qammod(0:M_data-1, M_data, 'UnitAveragePower', true)).^2) / kDataBits / Nfft;

interp_methods = {'linear', 'nearest', 'spline'};
numMethods     = length(interp_methods);
MSE_interp     = zeros(numMethods, length(fd_values), length(SNR_dB), length(profiles));

for idxP = 1:length(profiles)
    prof = profiles{idxP};
    fprintf('\nКанал %s (CP = %d)\n', prof.name, CP);

    for idxF = 1:length(fd_values)
        fd = fd_values(idxF);
        fprintf('  f_d = %d Гц (%d прогонов × %d символов)\n', ...
                fd, numRealizations, NumOfdmSymbols);
        
        % накопитель ошибок errBits(r, k) и MSE msePerRe(r, k)
        errBits  = zeros(numRealizations, length(SNR_dB));
        msePerRe = zeros(numRealizations, length(SNR_dB), numMethods);
        
        % внешний цикл - независимые прогоны
        for r = 1:numRealizations

            %% матрицыу
            bits_tx_data = randi([0 1], nDataSymbols * kDataBits, NumOfdmSymbols);
            mod_data = qammod(bits_tx_data, M_data, 'UnitAveragePower', true, 'InputType', 'bit');

            bits_tx_pil = randi([0 1], nPilots * kPilotBits, NumOfdmSymbols);
            mod_pilots  = qammod(bits_tx_pil, M_pilot, 'UnitAveragePower', true, 'InputType', 'bit');
            
            % OFDM-символы [Nfft x NumOfdmSymbols]
            Sp = zeros(Nfft, NumOfdmSymbols);
            Sp(all_subc(data_idx),  :) = mod_data;
            Sp(all_subc(pilot_idx), :) = mod_pilots;

            Sig = ifft(ifftshift(Sp, 1), Nfft, 1);
            
            % CP по столбцам
            Tx_Signal = [Sig(end-CP+1:end, :); Sig];

            %% канал  
            rayleighchan = comm.RayleighChannel(...
                'SampleRate',          Fs, ...
                'PathDelays',          prof.delays, ...
                'AveragePathGains',    prof.gains, ...
                'MaximumDopplerShift', fd, ...
                'NormalizePathGains',  true);

            Rx_vec = rayleighchan(Tx_Signal(:));
            Rx_Signal_clean = reshape(Rx_vec, size(Tx_Signal));
            release(rayleighchan);

            %% сигнал без шума, идеал H 
            Sig_clean = Rx_Signal_clean(CP+1 : CP+Nfft, :);
            Sp_clean  = fftshift(fft(Sig_clean, Nfft, 1), 1);
            H_true    = Sp_clean(all_subc, :) ./ Sp(all_subc, :);

            %% SNR
            for k = 1:length(SNR_dB)

                % шум
                N0 = Eb / (10^(SNR_dB(k)/10));
                noise = sqrt(N0/2) * (randn(size(Rx_Signal_clean)) + 1j*randn(size(Rx_Signal_clean)));
                Rx_Signal = Rx_Signal_clean + noise;
                
                % снятие CP + FFT по столбцам
                Sig_rec = Rx_Signal(CP+1 : CP+Nfft, :);
                Sp_rec  = fftshift(fft(Sig_rec, Nfft, 1), 1);

    
                % оценка по пилотам
                H_pilots = Sp_rec(all_subc(pilot_idx), :) ./ mod_pilots;
                H_pilots_mean = mean(H_pilots, 2);
                H_true_mean = mean(H_true, 2);

                % MSE цикл по инетрполяциям
                for m = 1:numMethods
                    H_est_m = interp1(pilot_idx(:), H_pilots_mean, (1:NumSC).', ...
                                      interp_methods{m}, 'extrap');
                    msePerRe(r, k, m) = channel_mse(H_est_m, H_true_mean, data_idx);
                end
                

                % ZF
                H_est = interp1(pilot_idx(:), H_pilots, 1:NumSC, 'linear', 'extrap');
                Sp_eq = Sp_rec(all_subc(data_idx), :) ./ H_est(data_idx, :);

                % демодуляция
                rx_bits = qamdemod(Sp_eq, M_data, 'UnitAveragePower', true, 'OutputType', 'bit');

                errBits(r, k) = sum(bits_tx_data(:) ~= rx_bits(:));
            end
        end

        BER_res(idxF, :, idxP) = sum(errBits, 1) / totalBits;

        idxSpline = find(strcmp(interp_methods, 'spline'));
        MSE_res(idxF, :, idxP) = mean(msePerRe(:, :, idxSpline), 1);

        for m = 1:numMethods
            MSE_interp(m, idxF, :, idxP) = mean(msePerRe(:, :, m), 1);
        end
    end
end

%% ДОП п.1 - ИХ и ЧХ каналов
for idxP = 1:length(profiles)
    prof = profiles{idxP};
    h_imp = zeros(1, round(max(prof.delays)*Fs) + 2);
    for l = 1:length(prof.delays)
        n = round(prof.delays(l) * Fs) + 1;
        h_imp(n) = h_imp(n) + 10^(prof.gains(l)/20);
    end

    idx = find(h_imp > 0);
    t_imp = (idx - 1) / Fs * 1e6;

    figure(1);
    subplot(3,2,(idxP-1)*2 + 1);
    stem(t_imp, h_imp(idx), 'filled', 'LineWidth', 1.5);
    grid on; xlabel('задержка, нс'); ylabel('амплитуда');
    title([prof.name ': ИХ (лин. масштаб)']); ylim([0 2]);

    subplot(3,2,(idxP-1)*2 + 2);
    stem(t_imp, 20*log10(h_imp(idx)), 'filled', 'LineWidth', 1.5);
    grid on; xlabel('задержка, нс'); ylabel('|h|, дБ');
    title([prof.name ': ИХ (в дБ)']); ylim([-40 5]);

    H = fftshift(fft(h_imp, Nfft));
    H = H / max(abs(H));
    f_axis = (-Nfft/2:Nfft/2-1) * (Fs/Nfft) / 1e6;

    figure(2); subplot(3,1,idxP);
    plot(f_axis, 20*log10(abs(H)+eps), 'LineWidth', 1.2);
    grid on; xlabel('f, МГц'); ylabel('|H|, дБ');
    title(sprintf('%s: АЧХ (fd = 0 Гц)', prof.name));
    xlim([-Fs/2e6 Fs/2e6]); ylim([-40 5]);
end

%% ДОП: созвездие
fd_demo   = 0;
snr_demo  = 20;

for idxP_demo = 1:length(profiles)
    prof_demo = profiles{idxP_demo};

    rayleighchan_demo = comm.RayleighChannel(...
        'SampleRate',          Fs, ...
        'PathDelays',          prof_demo.delays, ...
        'AveragePathGains',    prof_demo.gains, ...
        'MaximumDopplerShift', fd_demo, ...
        'NormalizePathGains',  true);

    bits_tx  = randi([0 1], nDataSymbols * kDataBits, 1);
    mod_data = qammod(bits_tx, M_data, 'UnitAveragePower', true, 'InputType', 'bit');

    bits_pilot  = randi([0 1], nPilots * kPilotBits, 1);
    mod_pilots  = qammod(bits_pilot, M_pilot, 'UnitAveragePower', true, 'InputType', 'bit');

    Sp = zeros(Nfft, 1);
    Sp(all_subc(data_idx))  = mod_data(:);
    Sp(all_subc(pilot_idx)) = mod_pilots(:);

    Sig = ifft(ifftshift(Sp), Nfft);
    Tx_Signal = [Sig(end-CP+1:end); Sig];

    Rx_Signal = rayleighchan_demo(Tx_Signal);
    release(rayleighchan_demo);

    N0_demo = Eb / (10^(snr_demo/10));
    noise   = sqrt(N0_demo/2) * (randn(size(Rx_Signal)) + 1j*randn(size(Rx_Signal)));
    Rx_Signal = Rx_Signal + noise;

    Sig_rec = Rx_Signal(CP+1 : CP+Nfft);
    Sp_rec  = fftshift(fft(Sig_rec, Nfft));

    H_pilots = Sp_rec(all_subc(pilot_idx)) ./ mod_pilots(:);
    H_est    = interp1(pilot_idx(:), H_pilots(:), (1:NumSC).', 'pchip', 'extrap');

    rx_before = Sp_rec(all_subc(data_idx));
    rx_after  = rx_before ./ H_est(data_idx);

    figure('Name', sprintf('созвездие %s', prof_demo.name), 'Color', 'w');

    subplot(1,2,1);
    plot(real(rx_before), imag(rx_before), 'b.', 'MarkerSize', 6);
    hold on;
    plot(real(mod_data), imag(mod_data), 'r+', 'MarkerSize', 8, 'LineWidth', 1);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title('до ZF');

    subplot(1,2,2);
    plot(real(rx_after), imag(rx_after), 'b.', 'MarkerSize', 6);
    hold on;
    plot(real(mod_data), imag(mod_data), 'r+', 'MarkerSize', 8, 'LineWidth', 1);
    grid on; axis equal;
    xlabel('I'); ylabel('Q');
    title('после ZF');

    sgtitle(sprintf('КАМ-%d: %s, f_d = %d Гц, SNR = %d дБ', ...
            M_data, prof_demo.name, fd_demo, snr_demo));
end

%% ДОП: п. 2 - BER для EPA/EVA/ETU fd=0
idxF_zero = find(fd_values == 0, 1);
figure(); hold on; grid on;
chan_markers = {'k-o', 'b-s', 'r-d'};
chan_colors  = {'k', 'b', 'r'};
for idxP = 1:length(profiles)
    prof = profiles{idxP};
    Y = reshape(BER_res(idxF_zero, :, idxP), [], 1);
    semilogy(SNR_dB, Y, chan_markers{idxP}, 'Color', chan_colors{idxP}, ...
             'LineWidth', 1.5, 'MarkerSize', 8, ...
             'DisplayName', sprintf('%s (f_d = 0 Гц)', prof.name));
end
xlabel('SNR (дБ)'); ylabel('BER');
title(sprintf('BER EPA/EVA/ETU при f_d = 0 Гц (КАМ-%d, CP = %d)', M_data, CP));
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-3 1]); hold off;

%% ДОП: п.3 - доплер спектры EPA5 EVA70 ETU300
fd_spec = [5, 70, 300];
names_spec = {'EPA', 'EVA', 'ETU'};
Fs_gains = Fs / (Nfft + CP);
figure();
for i = 1:length(names_spec)
    prof = profiles{i};
    fd = fd_spec(i);
    chan_dop = comm.RayleighChannel('SampleRate', Fs, ...
        'PathDelays', prof.delays, 'AveragePathGains', prof.gains, ...
        'MaximumDopplerShift', fd, 'NormalizePathGains', true, ...
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
    xlim([-1.5*fd, 1.5*fd]); ylim([-40, 5]);
    title(sprintf('%s, f_d = %d Гц', prof.name, fd));
    xlabel('f, Гц'); ylabel('норм. мощность, дБ'); grid on;
end

%% ДОП: п.4 - BER и Tc
markers = {'k-o', 'b-s', 'r-d', 'g-^', 'm-x'};
colors  = {'k', 'b', 'r', 'g', 'm'};
idxP = 1;
prof = profiles{idxP};
fd_Tc = [50, 2000];
T_ofdm_sym = (Nfft + CP) / Fs;

fprintf('\nTc для канала %s:\n', prof.name);
for fd = fd_Tc
    Tc = 1 / fd; Tc_ms = Tc * 1e3; Tc_sym = Tc / T_ofdm_sym;
    fprintf('  f_d = %4d Гц: Tc = %.3f мс (%.2f символов OFDM)\n', fd, Tc_ms, Tc_sym);
end

figure; hold on; grid on;
for idxF = 1:length(fd_values)
    fd = fd_values(idxF);
    Y  = reshape(BER_res(idxF, :, idxP), [], 1);
    semilogy(SNR_dB, Y, markers{idxF}, 'Color', colors{idxF}, ...
             'LineWidth', 1.5, 'MarkerSize', 8, ...
             'DisplayName', sprintf('f_d = %d Гц', fd));
end
xlabel('SNR (дБ)'); ylabel('BER');
title(sprintf('OFDM %s (CP = %d, КАМ-%d, %d прогонов на %d симв.)', ...
      prof.name, CP, M_data, numRealizations, NumOfdmSymbols));
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
ylim([1e-3 1]);
hold off;

%% ДОП: п.5 - MSE
figure(); hold on; grid on;
idxP_show = 1;
idxF_show = 1;
interp_markers = {'k-o', 'b-s', 'r-d', 'g-^', 'm-x'};

for m = 1:numMethods
    Y = reshape(MSE_interp(m, idxF_show, :, idxP_show), [], 1);
    semilogy(SNR_dB, Y, interp_markers{m}, 'LineWidth', 1.8, 'MarkerSize', 8, ...
             'DisplayName', interp_methods{m});
end

xlabel('SNR (дБ)'); ylabel('MSE');
title(sprintf('MSE оценки канала: %s, f_d = %d Гц (разные интерполяции)', ...
      profiles{idxP_show}.name, fd_values(idxF_show)));
legend('Location', 'southwest');
set(gca, 'YScale', 'log');
hold off;

%% MSE
function mse = channel_mse(H_est, H_true, idx)
    err = H_est(idx, :) - H_true(idx, :);
    mse = mean(abs(err(:)).^2);
end