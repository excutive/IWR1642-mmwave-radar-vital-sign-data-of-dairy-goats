function [process_adc, fs] = readDCA1000(Filename, process_num)
% readDCA1000 读取 DCA1000 采集的原始 bin 文件，去交错并取出 RX1 通道，慢时间抽稀
% 输入：
%   Filename    - bin 文件路径
%   process_num - 处理的帧数（默认 1024）
% 输出：
%   process_adc - numADCSamples x (numChirps/f_chirps) 的复数矩阵，RX1 通道抽稀后
%   fs          - 慢时间采样率（Hz）

%% 采集参数
numADCSamples = 256; % 每个 chirp 的 ADC 采样点数
numADCBits = 16;     % 每个采样点的 ADC 位数
numTX = 1;           % 发射天线数
numRX = 4;           % 接收天线数
numLanes = 2;        % LVDS 通道数，固定为 2，不要改
isReal = 0;          % 0 表示复采样数据，1 表示实采样
chirpLoop = 2;       % 每帧 chirp 数
f_chirps = 2;        % 抽稀因子：每 f_chirps 个 chirp 保留 1 个

if nargin < 2 || isempty(process_num)
    process_num = 1024; % 只对 process_num 帧的数据做处理
end

%% 读取 bin 文件
fid = fopen(Filename, 'r');
adcDataRow = fread(fid, 'int16');
if numADCBits ~= 16
    l_max = 2^(numADCBits-1)-1;
    adcDataRow(adcDataRow > l_max) = adcDataRow(adcDataRow > l_max) - 2^numADCBits;
end
fclose(fid);

fileSize = process_num*2*numADCSamples*numTX*numRX*2;
PRTnum = fix(fileSize/(numADCSamples*numRX));
fileSize = PRTnum * numADCSamples*numRX;
adcData = adcDataRow(1:fileSize);

%% 复采样（LVDS 去交错）
numChirps = fileSize/2/numADCSamples/numRX;
LVDS = zeros(1, fileSize/2);
counter = 1;
for i = 1:4:fileSize-1
    LVDS(1, counter)   = adcData(i)   + sqrt(-1)*adcData(i+2);
    LVDS(1, counter+1) = adcData(i+1) + sqrt(-1)*adcData(i+3);
    counter = counter + 2;
end

LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
LVDS = LVDS.';

%% 重组数据，取出 RX1
adcData = zeros(numRX, numChirps*numADCSamples);
for row = 1:numRX
    for i = 1:numChirps
        adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = ...
            LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
    end
end
retVal = reshape(adcData(1, :), numADCSamples, numChirps);

%% 慢时间抽稀：每 f_chirps 个 chirp 保留 1 个
process_adc = zeros(numADCSamples, numChirps/f_chirps);
for nchirp = 1:f_chirps:numChirps
    process_adc(:, (nchirp-1)/f_chirps+1) = retVal(:, nchirp);
end

% 慢时间采样率：帧周期 40 ms，每帧 chirpLoop 个 chirp，再抽稀 f_chirps
frame_period = 40e-3;
fs = (1/frame_period) * chirpLoop / f_chirps;  % 50 Hz / 2 = 25 Hz
end
