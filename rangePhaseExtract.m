function [angle_fft_last2, best_range_bin, deltaR] = rangePhaseExtract(process_adc)
% rangePhaseExtract 距离维 FFT、静态杂波抑制、最大能量选距离单元、相位提取与中心差分
% 输入：
%   process_adc - numADCSamples x numChirps 的 RX1 复数矩阵（readDCA1000 的输出）
% 输出：
%   angle_fft_last2 - 相位差分信号
%   best_range_bin   - 选中的距离单元
%   deltaR           - 距离分辨率（m）

%% 雷达参数
Fs = 10e6;          % ADC 采样率
c = 3e8;            % 光速
slope = 60.012e12;  % 调频斜率
startFreq = 77e9;   % 起始频率
numADCSamples = size(process_adc, 1);
RangFFT = 256;

%% 距离维 FFT
window = hamming(numADCSamples);
adcdata = process_adc;
numChirps = size(adcdata, 2);
for i = 1:numChirps
    adcdata(:, i) = adcdata(:, i) .* window;
end
fft_data = fft(adcdata, RangFFT);
fft_data = fft_data.';

%% 静态干扰消除
% 慢时间均值相减（去直流），前 5 个距离单元置零抑制近场泄漏
mean_fft = mean(fft_data, 1);
fft_data_filtered = fft_data - repmat(mean_fft, numChirps, 1);
fft_data_abs = abs(fft_data_filtered);
fft_data_abs(:, 1:5) = 0;

%% 最大能量准则选距离单元
deltaR = Fs*c/2/slope/RangFFT;
E_m = sum(fft_data_abs.^2, 1) / numChirps;  % 每个距离单元的累计能量
[~, best_range_bin] = max(E_m);

%% 相位提取
real_data = real(fft_data_filtered);
imag_data = imag(fft_data_filtered);
angle_fft_last = atan2(imag_data(:, best_range_bin), real_data(:, best_range_bin));
angle_fft_last = unwrap(angle_fft_last);

%% 中心差分
angle_fft_last2 = zeros(size(angle_fft_last));
for i = 2:numChirps-1
    angle_fft_last2(i) = (angle_fft_last(i+1) - angle_fft_last(i-1)) / 2;
end
angle_fft_last2(1) = angle_fft_last(2) - angle_fft_last(1);
angle_fft_last2(end) = angle_fft_last(end) - angle_fft_last(end-1);
end
