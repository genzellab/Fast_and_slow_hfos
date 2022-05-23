function [S, E, M] = findHFOs(Filt_EEG, timestamps, DetectThreshold, LimitThreshold,varargin)
% [S, E, M] = findRipples(Filt_EEG, DetectThreshold, LimitThreshold, varargin)
% 
% INPUTS: 
% Filt_EEG: EEG tsd filtered in the ripples (e.g. 100-300 Hz) range
% DetectThreshold: Threshold for detection of ripples 
% LimitThreshold: Threshold for finding the ripple boundaries
%
% PARAMETERS:
% Q1: number of cycles to check to find boundaries
% CloseThreshold: Closeness threhshold,if two ripple events are closer than
%   this,  lump them together
% MinRippleDuration: discard theta evants shorter than this
% OUTPUTS:
% S: ripple events start times
% E: ripple events end times
% M: ts object with ripple events peak times

% copyright (c) 1999 Francesco P. Battaglia
% This software is released under the GNU GPL
% www.gnu.org/copyleft/gpl.html


% adapted for LG data FPB 2016
% adapted for LG data AAZ 2020

% parameters
Q1 = 3;
%CloseThreshold = 20 / 1000;
CloseThreshold = 200 / 1000;

MinRippleDuration = 30 / 1000;
MaxRippleDuration = 100 / 1000;


EEGStart = timestamps(1);
EEGEnd = timestamps(end);

% do it in chunks
% chunk length in seconds

LChunk = 300; %Might have to change

ChStart = EEGStart;
ChEnd = min(ChStart + LChunk, EEGEnd); %ChStart should thus be given in seconds. 
nChunks = (EEGEnd - EEGStart)/LChunk;
TRStart = [];
TREnd = [];
TRMax = [];
TRMaxValue = [];

i = 0;
% % % % % % % % h = waitbar(0, 'Find Ripples...');

while ChStart < EEGEnd
% % % % % % % % % %   waitbar(i/nChunks, h);

vecount=[1:length(timestamps)].';

t1=vecount(timestamps >= ChStart);
t1=t1(1); 

t2=vecount(timestamps >= ChEnd);
t2=t2(1); 

%   t1 = find(timestamps >= ChStart, 1);
%   t2 = find(timestamps >= ChEnd, 1);
  t = timestamps(t1:t2);
  eeg = Filt_EEG(t1:t2);
  
  sz = size(eeg);
  if sz(1) ~= 1
    eeg = eeg';    % we want to work with row vectors
    t = t';
  end
  
  

  
  de = diff(eeg);
  de1 = [de 0];
  de2 = [0 de];
  
  
  %finding peaks

  vecount=vecount.';
  
  upPeaksIdx=vecount(de1 < 0 & de2 > 0);
  downPeaksIdx = vecount(de1 > 0 & de2 < 0);
  
%   upPeaksIdx = find(de1 < 0 & de2 > 0);
%   downPeaksIdx = find(de1 > 0 & de2 < 0);
  
  PeaksIdx = [upPeaksIdx downPeaksIdx];
  PeaksIdx = sort(PeaksIdx);
  
  Peaks = eeg(PeaksIdx);
  Peaks = abs(Peaks);
  
  %when a peaks is above threshold, we detect a ripple
    vecount2=[1:length(Peaks)];
  RippleDetectIdx=vecount2(Peaks > DetectThreshold);
%   RippleDetectIdx = find(Peaks > DetectThreshold);

  DetectDiff = [0 diff(RippleDetectIdx)];
  RippleDetectIdx = RippleDetectIdx(DetectDiff > 2);
  RippleDetectIdx = RippleDetectIdx(RippleDetectIdx < length(Peaks)-Q1+1);
  RippleStart = zeros(1, length(RippleDetectIdx));
  RippleEnd = zeros(1, length(RippleDetectIdx));
  RippleMax = zeros(1, length(RippleDetectIdx));
  RippleMaxValue = zeros(1, length(RippleDetectIdx));
  for ii = 1:length(RippleDetectIdx) 
    CP = RippleDetectIdx(ii); % Current Peak
    % detect start of the ripple
    for j = CP-1:-1:Q1
      if all(Peaks(j-Q1+1:j) < LimitThreshold)
	
	break;
      end
    end
    RippleStart(ii) = j;
    %detect end of ripple
    for j = CP+1:length(Peaks)-Q1+1
      if all(Peaks(j:j+Q1-1)< LimitThreshold)
	
	break;
      end
    end
    RippleEnd(ii) = j;
    [RippleMaxValue(ii), RippleMax(ii)] = max(Peaks(RippleStart(ii):RippleEnd(ii)));
    RippleMax(ii) = RippleStart(ii) + RippleMax(ii) - 1;
  end
  
  TRStart = [TRStart t(PeaksIdx(RippleStart))]; %#ok<*AGROW>
  TREnd = [TREnd t(PeaksIdx(RippleEnd))];
  TRMax = [TRMax t(PeaksIdx(RippleMax))];
  TRMaxValue = [TRMaxValue RippleMaxValue];
  i = i+1;
  ChStart = ChStart + LChunk;

  ChEnd = min(ChStart + LChunk, EEGEnd);
end
% % % % % % % % % close(h);
i = 2;


while i <= length(TRStart)
  
  if (TRStart(i) - TREnd(i-1)) < CloseThreshold
    TRStart = [TRStart(1:i-1) TRStart(i+1:end)];
    TREnd = [TREnd(1:i-2) TREnd(i:end)];
    [~, ix] = max([TRMaxValue(i-1) TRMaxValue(i)]);
    TRMax = [TRMax(1:i-2) TRMax(i - 2 + ix) TRMax(i+1:end)]; 
  else
    i = i +1 ;
    
  end
  
end
length(TRStart);
length(TREnd);
i = 1;
while i <= length(TRStart)
  if(TREnd(i) - TRStart(i)) < MinRippleDuration
    TRStart = [TRStart(1:i-1) TRStart(i+1:end)];
    TREnd = [TREnd(1:i-1) TREnd(i+1:end)];
    TRMax = [TRMax(1:i-1), TRMax(i+1:end)];
  else
    i = i+ 1;
    
  end
  
end

%%Max duration criteria
i = 1;
while i <= length(TRStart)
  if(TREnd(i) - TRStart(i)) >= MaxRippleDuration
    TRStart = [TRStart(1:i-1) TRStart(i+1:end)];
    TREnd = [TREnd(1:i-1) TREnd(i+1:end)];
    TRMax = [TRMax(1:i-1), TRMax(i+1:end)];
  else
    i = i+ 1;
    
  end
  
end



S = (TRStart);
E = (TREnd);
M = (TRMax);