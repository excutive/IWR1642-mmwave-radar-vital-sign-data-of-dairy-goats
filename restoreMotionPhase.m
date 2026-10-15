function angle_fft_clean = restoreMotionPhase(angle_fft_last2, motion_mask, fs, min_valid)
% restoreMotionPhase 体动恢复：频-相轨迹重建
% 对每个体动区间，用左右相邻有效段的主频线性插值瞬时频率，再积分重建相位差分
% 输入：
%   angle_fft_last2 - 相位差分信号（列向量）
%   motion_mask     - 体动掩码（true 表示体动）
%   fs              - 慢时间采样率（Hz）
%   min_valid       - 有效段的最短长度（样本数）
% 输出：
%   angle_fft_clean - 恢复后的相位差分信号

if nargin < 4 || isempty(min_valid), min_valid = 16; end

angle_fft_last2 = angle_fft_last2(:);
motion_mask = motion_mask(:);
angle_fft_clean = angle_fft_last2;

%% 找体动区间
d = diff([false; motion_mask; false]);
s_list = find(d == 1);
e_list = find(d == -1) - 1;

%% 全局中位频率：所有足够长的有效段主频的中位数
valid = ~motion_mask;
d = diff([false; valid; false]);
vs = find(d == 1);
ve = find(d == -1) - 1;
seg_freqs = [];
for k = 1:numel(vs)
    if ve(k) - vs(k) + 1 >= min_valid
        seg_freqs(end+1) = dominantFreq(angle_fft_last2(vs(k):ve(k)), fs); %#ok<AGROW>
    end
end
if isempty(seg_freqs)
    f_g = 0.25;
else
    f_g = median(seg_freqs);
end

%% 逐个体动区间重建
for k = 1:numel(s_list)
    s = s_list(k);
    e = e_list(k);
    % 左右相邻有效段的主频，缺失时用全局中位频率
    if s - 1 >= min_valid
        f_L = dominantFreq(angle_fft_last2(s-min_valid:s-1), fs);
    else
        f_L = f_g;
    end
    if numel(angle_fft_last2) - e >= min_valid
        f_R = dominantFreq(angle_fft_last2(e+1:e+min_valid), fs);
    else
        f_R = f_g;
    end
    T = e - s;
    for t = s:e
        if T > 0
            ft = f_L + (f_R - f_L) * (t - s) / T;   % 线性插值瞬时频率
        else
            ft = f_L;
        end
        angle_fft_clean(t) = 2*pi*ft/fs;            % 相位差分 ≈ 2πf/fs
    end
end
end


function f0 = dominantFreq(x, fs)
% dominantFreq 求信号在 [0.1, 2.0] Hz 内的 FFT 主频
N = numel(x);
spec = abs(fft(x));
f_axis = (0:N-1) * (fs/N);
band = (f_axis >= 0.1) & (f_axis <= 2.0);
f_band = f_axis(band);
spec_band = spec(band);
if isempty(f_band)
    f0 = 0.25;
else
    [~, idx] = max(spec_band);
    f0 = f_band(idx);
end
end
