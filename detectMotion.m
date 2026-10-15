function motion_mask = detectMotion(angle_fft_last2, fs, Q_hi, Q_lo, w, stftWin, stftHop, movWin, eps_std)
% detectMotion 体动检测：7 谱特征 + 稳健标准化 + 双阈值迟滞
% 输入：
%   angle_fft_last2 - 相位差分信号（列向量）
%   fs              - 慢时间采样率（Hz）
%   Q_hi, Q_lo      - 迟滞系数
%   w               - 7 维特征权重
%   stftWin, stftHop- STFT 帧长/步长（可调）
%   movWin          - 稳健标准化的滑动中位数窗口（可调）
%   eps_std         - 标准化除零小量
% 输出：
%   motion_mask     - 与输入等长的逻辑列向量，true 表示体动

%% 默认参数
if nargin < 3 || isempty(Q_hi), Q_hi = 3; end
if nargin < 4 || isempty(Q_lo), Q_lo = 1.5; end
if nargin < 5 || isempty(w), w = [1 1 1 1 -1 1 1]; end  % 峰度随运动下降，取负权
if nargin < 6 || isempty(stftWin), stftWin = 64; end
if nargin < 7 || isempty(stftHop), stftHop = 8; end
if nargin < 8 || isempty(movWin), movWin = 9; end
if nargin < 9 || isempty(eps_std), eps_std = 0.01; end

angle_fft_last2 = angle_fft_last2(:);

%% STFT，功率谱 P = |S|^2
[S, f] = stft(angle_fft_last2, fs, 'Window', hann(stftWin, 'periodic'), ...
    'OverlapLength', stftWin - stftHop, 'FFTLength', stftWin);
P = abs(S).^2;
n_frames = size(P, 2);

%% 逐帧计算 7 个谱形态特征
feat = zeros(7, n_frames);
for j = 1:n_frames
    Pj = P(:, j) + 1e-12;
    p = Pj / sum(Pj);
    feat(1, j) = exp(mean(log(Pj))) / mean(Pj);              % 谱平坦度
    feat(2, j) = -sum(p(p > 0) .* log(p(p > 0)));            % 谱熵
    centroid = sum(f .* p);                                  % 谱质心
    feat(4, j) = sqrt(sum((f - centroid).^2 .* p));          % 谱带宽
    feat(5, j) = max(Pj) / mean(Pj);                         % 谱峰度 crest
    cum = cumsum(Pj);                                        % 95% 滚降频率
    idx = find(cum >= 0.95*cum(end), 1);
    feat(6, j) = f(idx);
end
% 谱通量、delta log 能量需要相邻帧
for j = 2:n_frames
    feat(3, j) = sum((abs(S(:, j)) - abs(S(:, j-1))).^2);
    feat(7, j) = abs(log(sum(P(:, j)) + 1e-12) - log(sum(P(:, j-1)) + 1e-12));
end

%% 稳健标准化
z = zeros(size(feat));
for i = 1:7
    med = movmedian(feat(i, :), movWin);
    mad_v = movmad(feat(i, :), movWin);
    z(i, :) = (feat(i, :) - med) ./ (mad_v + eps_std);
end

%% 综合判别得分
comp_score = w * z;

%% 双阈值迟滞
c_med = median(comp_score);
c_mad = median(abs(comp_score - c_med));   % 中位数绝对偏差（MAD）
thr_hi = c_med + Q_hi * c_mad;
thr_lo = c_med + Q_lo * c_mad;

strong = comp_score > thr_hi;
weak = comp_score > thr_lo;
frame_motion = strong;
% 迟滞扩展：与已检出的运动帧相邻、且大于 thr_lo 的帧并入
for iter = 1:n_frames
    dilated = [false, frame_motion(1:end-1)] | [frame_motion(2:end), false];
    new = dilated & weak & ~frame_motion;
    if ~any(new), break; end
    frame_motion = frame_motion | new;
end

%% 后处理：去孤立尖峰、填短空洞
frame_motion = fillGapsRemoveSpikes(frame_motion, 2, 1);

%% 帧级掩码映射回样本级
motion_mask = false(size(angle_fft_last2));
for j = find(frame_motion)
    s = (j-1)*stftHop + 1;
    e = min(s + stftWin - 1, numel(motion_mask));
    motion_mask(s:e) = true;
end
end


function m = fillGapsRemoveSpikes(m, gapFill, minRun)
% fillGapsRemoveSpikes 后处理：填短空洞、去孤立尖峰
% 去孤立尖峰：连续 true 长度小于 minRun 的置 false
d = diff([false, m, false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
for k = 1:numel(starts)
    if ends(k) - starts(k) + 1 < minRun
        m(starts(k):ends(k)) = false;
    end
end
% 填短空洞：连续 false 长度不超过 gapFill 且两侧为 true 的置 true
d = diff([false, ~m, false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
for k = 1:numel(starts)
    if ends(k) - starts(k) + 1 <= gapFill && starts(k) > 1 && ends(k) < numel(m)
        m(starts(k):ends(k)) = true;
    end
end
end
