 clc;
clear;
close all;
%%
FLAG_plot_constellations = false; 
FLAG_awgn_only  = true;

FLAG_cp_sync = true;   
FLAG_minn_sync = false; 

FLAG_time_shift = true;

FLAG_cfo = false;
FLAG_cfo_comp = false;   
cfo_list = [-0.5, - 0.3, -0.2, 0, 0.2, 0.3, 0.5];     

if FLAG_cp_sync && FLAG_minn_sync
    error('метод синхронизации');
end

Nfft = 2048;
T = 100e-6;
Fs = Nfft / T;

NumSC = 1600; 
NGI = Nfft - NumSC; 

GI_left = 1 : (NGI/2);               
GI_right = (Nfft - NGI/2 + 1) : Nfft; 
zeros_c = Nfft/2 + 1;  
exc_sc = [GI_left, GI_right, zeros_c];

all_sc = -Nfft/2 : (Nfft/2 - 1);
all_sc_shifted = all_sc + Nfft/2 + 1; 
int_sc = all_sc_shifted(~ismember(all_sc_shifted, exc_sc));
NumSC_actual = length(int_sc); 

M_data  = 16;
M_pilot = 4;
SNR_dB  = 0:2:20; 

NumOfdmSymbols  = 30; 
numRealizations = 1; 

pilot_step = 8;
kDataBits = log2(M_data);
kPilotBits = log2(M_pilot);
CP = Nfft/4;

tau_list = [0, -CP/2, -2*CP, CP/2, 2*CP];
tau_names = {'dt=0', ...
             sprintf('dt=-%d (|dt|<CP)', CP/2), ...
             sprintf('dt=-%d (|dt|>CP)', 2*CP), ...
             sprintf('dt=+%d (dt<CP)', CP/2), ...
             sprintf('dt=+%d (dt>CP)', 2*CP)};

fd_values = [0]; 
profiles = {
    struct('name','EPA', 'delays',[0 30 70 90 110 190 410]*1e-9, ...
           'gains',[0 -1 -2 -3 -8 -17.2 -20.8]);
    struct('name','EVA', 'delays',[0 30 150 310 370 710 1090 1730 2510]*1e-9, ...
           'gains',[0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9]);
    struct('name','ETU', 'delays',[0 50 120 200 230 500 1600 2300 5000]*1e-9, ...
           'gains',[-1 -1 -1 0 0 0 -3 -5 -7]);
};

pilot_idx = 1:pilot_step:NumSC_actual;
data_idx = setdiff(1:NumSC_actual, pilot_idx);
nPilots = length(pilot_idx);
nDataSymbols = length(data_idx);

totalBits = numRealizations * (NumOfdmSymbols - 1) * nDataSymbols * kDataBits;
Eb = mean(abs(qammod(0:M_data-1, M_data, 'UnitAveragePower', true)).^2) / kDataBits / Nfft;

metric_saved = [];   
delta_saved  = 0;

%% преамбула
L_minn = Nfft / 4;
A_minn = qammod(randi([0 3], L_minn, 1), 4, 'UnitAveragePower', true); 
Minn_Preamble = [A_minn; A_minn; -A_minn; -A_minn];   

%%
if ~FLAG_time_shift
    tau_list_use = [0];
    tau_names_use = {'dt=0'};
else

    if FLAG_minn_sync
        tau_list_use = [0, CP/2, 2*CP];
        tau_names_use = {'dt=0', ...
                         sprintf('dt=+%d (dt<CP)', CP/2), ...
                         sprintf('dt=+%d (dt>CP)', 2*CP)};
    else
        tau_list_use = tau_list;
        tau_names_use = tau_names;
    end
end

if FLAG_minn_sync
    sync_name = 'Minn';
    preamble_skip = Nfft;
else
    sync_name     = 'CP';
    preamble_skip = 0;
end

if FLAG_awgn_only
    num_profiles_to_run = 1; 
else
    num_profiles_to_run = length(profiles);
end

if FLAG_cfo
    cfo_list_use = cfo_list;
else
    cfo_list_use = 0;
end

for idxP = 1:num_profiles_to_run
    prof = profiles{idxP};
    if FLAG_awgn_only
        prof.name = 'AWGN';
    end
    fprintf('%s ( %s)\n', prof.name, sync_name);

    BER_res = zeros(length(tau_list_use), length(SNR_dB), length(cfo_list_use));
    BER_res_noeq = zeros(length(tau_list_use), length(SNR_dB), length(cfo_list_use));

    for idxCFO = 1:length(cfo_list_use)
        cfo_norm = cfo_list_use(idxCFO);
        for idxCase = 1:length(tau_list_use)
            delta_true = tau_list_use(idxCase);
            fprintf('%s, CFO=%.2f\n', tau_names_use{idxCase}, cfo_norm);

            errBits = zeros(numRealizations, length(SNR_dB));
            errBits_noeq = zeros(numRealizations, length(SNR_dB));

            for k = 1:length(SNR_dB)
                fprintf('SNR %d дБ\n', SNR_dB(k));
                for r = 1:numRealizations
                    %% данные
                    bits_tx_data = randi([0 1], nDataSymbols*kDataBits, NumOfdmSymbols);
                    mod_data = qammod(bits_tx_data, M_data, 'UnitAveragePower', true, 'InputType', 'bit');

                    bits_tx_pil = randi([0 1], nPilots*kPilotBits, NumOfdmSymbols);
                    mod_pilots = qammod(bits_tx_pil, M_pilot, 'UnitAveragePower', true, 'InputType', 'bit');

                    Sp = zeros(Nfft, NumOfdmSymbols);
                    Sp(int_sc(data_idx),  :) = mod_data;
                    Sp(int_sc(pilot_idx), :) = mod_pilots;
                    Sig = ifft(ifftshift(Sp, 1), Nfft, 1);
                    Sig_with_cp = [Sig(end-CP+1:end, :); Sig];

                    %% TX
                    if FLAG_minn_sync
                        Tx_vec = [Minn_Preamble; Sig_with_cp(:)];
                    else
                        Tx_vec = Sig_with_cp(:);
                    end

                    %% канал
                    if FLAG_awgn_only
                        Rx_vec = Tx_vec;
                    else
                        rayleighchan = comm.RayleighChannel('SampleRate', Fs, ...
                            'PathDelays', prof.delays, ...
                            'AveragePathGains', prof.gains, ...
                            'MaximumDopplerShift', fd_values(1),...
                            'NormalizePathGains', true);

                        Rx_vec = rayleighchan(Tx_vec);
                        release(rayleighchan);
                    end

                    %% dt
                    if FLAG_time_shift && delta_true ~= 0
                        if delta_true > 0
                            Rx_vec = [zeros(delta_true, 1); Rx_vec(1:end-delta_true)];
                        elseif delta_true < 0
                            d = abs(delta_true);
                            Rx_vec = [Rx_vec(d+1:end); zeros(d, 1)];
                        end
                    end

                    %% df
                        if FLAG_cfo && cfo_norm ~= 0
                             n = (0:numel(Rx_vec)-1).';
                             Rx_vec = Rx_vec .* exp(1j*2*pi*cfo_norm*n/Nfft);
                        end

                    %% шум
                    N0 = Eb / (10^(SNR_dB(k)/10));
                    noise = sqrt(N0/2) * (randn(size(Rx_vec)) + 1j*randn(size(Rx_vec)));
                    Rx_noisy = Rx_vec + noise;

                    %% синхронизация
                    if FLAG_cp_sync
                        [delay_est, cfo_norm_est, ACF_cp] = estimate_timing_cp(Rx_noisy, Nfft, CP, NumOfdmSymbols-1);
                        
                    elseif FLAG_minn_sync
                        [delay_est, metric, cfo_norm_est] = estimate_timing_minn(Rx_noisy, Nfft, CP);
                        if idxP == 1 && idxCFO == 1 && idxCase == min(2, length(tau_list_use)) && k == length(SNR_dB) && r == 1
                            metric_saved = metric;
                            delta_saved  = delta_true;
                        end
                    else
                        delay_est = 0;
                        cfo_norm_est = 0;
                    end

                    %% компенсация cfo
                    if FLAG_cfo_comp && cfo_norm_est ~= 0
                        n = (0:numel(Rx_noisy)-1).';
                        Rx_noisy = Rx_noisy .* exp(-1j*2*pi*cfo_norm_est*n/Nfft);
                    end

                    % выравнивание (скомпенсированного сигнала)
                    rx_data_part = Rx_noisy(delay_est + preamble_skip + 1 : end);
                    rx_aligned = [rx_data_part; zeros(delay_est, 1)];
                    Rx_Matrix = reshape(rx_aligned, Nfft+CP, NumOfdmSymbols);

                    %% снятие CP + FFT (- послед. символ)
                    Sig_rec = Rx_Matrix(CP+1:CP+Nfft, 1:end-1);
                    Sp_rec = fftshift(fft(Sig_rec, Nfft, 1), 1);

                    %% эквалайзер + демодуляция
                    H_pilots = Sp_rec(int_sc(pilot_idx), :) ./ mod_pilots(:, 1:end-1);
                    H_est = interp1(pilot_idx(:), H_pilots, 1:NumSC_actual, 'linear', 'extrap');
                    Sp_eq = Sp_rec(int_sc(data_idx), :) ./ H_est(data_idx, :);

                    rx_bits = qamdemod(Sp_eq, M_data, 'UnitAveragePower', true, 'OutputType', 'bit');
                    rx_bits_noeq = qamdemod(Sp_rec(int_sc(data_idx), :), M_data, 'UnitAveragePower', true, 'OutputType', 'bit');

                    bits_tx_data_ber = bits_tx_data(:, 1:end-1);
                    errBits(r, k) = sum(bits_tx_data_ber(:) ~= rx_bits(:));
                    errBits_noeq(r, k) = sum(bits_tx_data_ber(:) ~= rx_bits_noeq(:));
                        %% созвездия
                        if FLAG_plot_constellations && idxCase == min(4, length(tau_list_use)) && k == length(SNR_dB) && r == 1
    
                           rx_before = Sp_rec(int_sc(data_idx), :);
                           scatterplot(rx_before(:)); 
                           title(sprintf('до ZF (%s, CFO=%.2f, dt=%d)', prof.name, cfo_norm, delta_true));
    
                            scatterplot(Sp_eq(:));
                            title(sprintf('после ZF (%s, CFO=%.2f, dt=%d)', prof.name, cfo_norm, delta_true));
                            drawnow;
                        end
                end
            end

            BER_res(idxCase, :, idxCFO) = sum(errBits, 1) / totalBits;
            BER_res_noeq(idxCase, :, idxCFO) = sum(errBits_noeq, 1) / totalBits;
        end
    end

    %% BER
    figure('Name', sprintf('BER (%s, %s)', ...
        prof.name, sync_name), 'Color', 'w');

    nCurves = length(tau_list_use) * length(cfo_list_use);

    colors = hsv(nCurves);

    markers = {'o', 's', '^', 'd', 'p', 'v', 'x', ...
               '+', '*', '<', '>', 'h', '.'};

    lineStyles = {'-', '--', '-.', ':'};

    for eqCase = 1:2

        subplot(1,2,eqCase);
        hold on;
        grid on;

        curveIdx = 0;

        for idxCFO = 1:length(cfo_list_use)
            for idxCase = 1:length(tau_list_use)

                curveIdx = curveIdx + 1;

                if FLAG_time_shift && FLAG_cfo
                    legend_name = sprintf('%s, df=%.2f', ...
                        tau_names_use{idxCase}, ...
                        cfo_list_use(idxCFO));

                elseif FLAG_cfo
                    legend_name = sprintf('df=%.2f', ...
                        cfo_list_use(idxCFO));

                elseif FLAG_time_shift
                    legend_name = tau_names_use{idxCase};

                else
                    legend_name = 'dt=0, df=0';
                end

                if eqCase == 1
                    ber = BER_res(idxCase, :, idxCFO);
                else
                    ber = BER_res_noeq(idxCase, :, idxCFO);
                end

                markerIdx = mod(curveIdx-1, ...
                    length(markers)) + 1;

                styleIdx = mod(floor((curveIdx-1) / ...
                    length(markers)), length(lineStyles)) + 1;

                semilogy(SNR_dB, ber, ...
                    'LineStyle', lineStyles{styleIdx}, ...
                    'Marker', markers{markerIdx}, ...
                    'Color', colors(curveIdx,:), ...
                    'LineWidth', 1.5, ...
                    'MarkerSize', 7, ...
                    'DisplayName', legend_name);

            end
        end

        xlabel('SNR (дБ)');
        ylabel('BER');

        if eqCase == 1
            title(sprintf('С эквалайзером (%s, %s)', ...
                prof.name, sync_name));
        else
            title(sprintf('Без эквалайзера (%s, %s)', ...
                prof.name, sync_name));
        end

        legend('Location', 'southwest');
        set(gca, 'YScale', 'log');
        ylim([1e-4 1]);
        hold off;

    end
end

%% метрика Минна
if FLAG_minn_sync && ~isempty(metric_saved)
    figure('Name', 'метрика Минна', 'Color', 'w');
    hold on; 
    grid on;
    plot(metric_saved, 'b', 'LineWidth', 1.2, 'DisplayName', 'Метрика');

    [~, max_idx] = max(metric_saved);
    plot(max_idx, metric_saved(max_idx), 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('%d', max_idx-1));

    xlabel('отсчеты'); 
    ylabel('метрика Минна M(d)');
    title(sprintf('АКФ Минна (SNR=%d дБ, сдвиг=%d, CFO=%.2f)', SNR_dB(end), delta_saved, cfo_list_use(1)));
    legend('Location', 'best');
    xlim([1, 3*Nfft]);
end

fprintf('\nвсё\n');

%% функции
function [delay_est, cfo_norm_est, ACF_total] = estimate_timing_cp(Rx_vec, Nfft, CP, NumSymsAvg)
    Rx_vec = Rx_vec(:);
    sym_len = Nfft + CP;
    ACF_total = zeros(Nfft, 1);
    max_syms = min(NumSymsAvg, floor(length(Rx_vec)/sym_len) - 1);

    for sym_idx = 1:max_syms
        start_idx = (sym_idx-1)*sym_len + 1;
        if (start_idx + 2*Nfft + CP - 2) <= length(Rx_vec)
            left_all  = Rx_vec(start_idx : start_idx + Nfft + CP - 2); 
            right_all = Rx_vec(start_idx + Nfft : start_idx + 2*Nfft + CP - 2); 
            prod = conj(left_all) .* right_all;
            ACF_sym = conv(prod, ones(CP, 1), 'valid');
            ACF_total = ACF_total + ACF_sym;
        end
    end
    [~, max_idx] = max(abs(ACF_total));
    delay_est = max_idx - 1; %кросс-корреляция
    
    %дельта f = (1/2pi) * arg(CF_max)
    cfo_norm_est = angle(ACF_total(max_idx)) / (2*pi); 
end


    % function [delay_est, metric, cfo_norm_est] = estimate_timing_minn(Rx_vec, Nfft, CP)
    %     Rx_vec = Rx_vec(:);
    %     L = Nfft / 4;
    %     search_len = length(Rx_vec) - 4*L;
    %     metric = zeros(search_len, 1);
    % 
    %     best_metric = -inf;
    %     P_at_max = 0;
    % 
    %     for d = 1:search_len
    %         left1  = Rx_vec(d : d + L - 1);
    %         right1 = Rx_vec(d + L : d + 2*L - 1);
    %         left2  = Rx_vec(d + 2*L : d + 3*L - 1);
    %         right2 = Rx_vec(d + 3*L : d + 4*L - 1);
    % 
    %         P = sum(conj(left1) .* right1) + sum(conj(left2) .* right2);
    %         R = sum(abs(left1).^2);
    % 
    %         if R > 0
    %             metric(d) = (abs(P)^2) / (R^2);
    %         end
    %         
    %         if metric(d) > best_metric
    %             best_metric = metric(d);
    %             P_at_max = P;
    %         end
    %     end
    % 
    %     [~, max_idx] = max(metric);
    %     delay_est = max_idx - 1; 
    %     
    %     cfo_norm_est = angle(P_at_max) * 2 / pi;
    % end

function [delay_est, metric, cfo_norm_est] = estimate_timing_minn(Rx_vec, Nfft, CP)
    Rx_vec = Rx_vec(:);
    L = Nfft / 4;
    search_len = length(Rx_vec) - 4*L;

    if search_len <= 0
        delay_est = 0; 
        metric = []; 
        cfo_norm_est = 0; 
        return;
    end
    %1 и 2
    v1 = conj(Rx_vec(1 : end-L)) .* Rx_vec(L+1 : end);
    %3 и 4
    v2 = conj(Rx_vec(2*L+1 : end-L)) .* Rx_vec(3*L+1 : end);

    % P1 P2
    P1_all = conv(v1, ones(L, 1), 'valid');
    P2_all = conv(v2, ones(L, 1), 'valid');

    % P(d)
    P = P1_all(1:search_len) + P2_all(1:search_len);

    metric = abs(P); 

    % пик
    [~, max_idx] = max(metric);
    delay_est = max_idx - 1; 
    
    % df
    cfo_norm_est = angle(P(max_idx)) * 2 / pi;
end