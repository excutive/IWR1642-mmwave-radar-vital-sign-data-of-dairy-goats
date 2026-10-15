function [f_resp, RR, best_idx] = aimfEmd(breath_data, fs, f_ref, eps_E, eps_f, eps_p, f_lo, f_hi)
% aimfEmd AIMF-EMD：EMD 分解 + 三判据自动选 IMF + periodogram 求呼吸率
% 输入：
%   breath_data - 自适应滤波后的呼吸信号（列向量）
%   fs          - 慢时间采样率（Hz）
%   f_ref       - 参考频率（自适应滤波阶段估计的主频）
%   eps_E       - 能量一致性的除零小量
%   eps_f       - 谱对齐的除零小量
%   eps_p       - 周期性的除零小量
%   f_lo, f_hi  - 生理频带（Hz）
% 输出：
%   f_resp      - 呼吸频率（Hz）
%   RR          - 呼吸率（次/分钟）
%   best_idx    - 选中的 IMF 序号

if nargin < 3 || isempty(f_ref), f_ref = 0.25; end
if nargin < 4 || isempty(eps_E), eps_E = 0.01; end
if nargin < 5 || isempty(eps_f), eps_f = 0.001; end
if nargin < 6 || isempty(eps_p), eps_p = 1e-3; end
if nargin < 7 || isempty(f_lo), f_lo = 0.1; end
if nargin < 8 || isempty(f_hi), f_hi = 2.0; end

breath_data = breath_data(:);
N = numel(breath_data);

%% EMD 分解
imf = emd(breath_data);
num_imfs = size(imf, 2);

%% 各 IMF 的能量
energies = sum(imf.^2, 1);
E_med = median(energies);
f_peak = zeros(1, num_imfs);
C_p = zeros(1, num_imfs);
f_axis = (0:N-1) * (fs/N);
band = (f_axis >= f_lo) & (f_axis <= f_hi);
f_band = f_axis(band);
for i = 1:num_imfs
    % 主频（Eq.35，在生理频带内）
    spec = abs(fft(imf(:, i)));
    [~, idx] = max(spec(band));
    f_peak(i) = f_band(idx);
    % 自相关峰值间隔的周期性（Eq.37-38）
    acf = xcorr(imf(:, i), 'coeff');
    acf = acf(N:end);
    [~, locs] = findpeaks(acf);
    if numel(locs) > 1
        C_p(i) = 1 / (std(diff(locs)) + eps_p);
    end
end

%% 三判据综合得分（Eq.34、Eq.36、Eq.39）
score = zeros(1, num_imfs);
for i = 1:num_imfs
    C_E = 1 / (abs(energies(i) - E_med) + eps_E);
    C_f = 1 / (abs(f_peak(i) - f_ref) + eps_f);
    score(i) = C_E * C_f * C_p(i);
end

%% 选呼吸分量（Eq.40）
[~, best_idx] = max(score);
IMF_resp = imf(:, best_idx);

%% periodogram 求最终频率（Eq.41-42）
[P, f] = periodogram(IMF_resp, hamming(N), [], fs);
band = (f >= f_lo) & (f <= f_hi);
f_band = f(band);
P_band = P(band);
[~, idx] = max(P_band);
f_resp = f_band(idx);
RR = f_resp * 60;
end
