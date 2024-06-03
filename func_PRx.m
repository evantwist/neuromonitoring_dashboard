%% Function for PRx analysis
function [PRx] = func_PRx(ICP, MAP)
% Input: 
%       ICP: Continuous intracranial pressure time series
%       MAP: Continuous mean arterial pressure time series
% Time series should be synchronized and span the same timeperiod, and
% optionally filtered to remove noise from pulse rate and respiration.
%
% Output: 
%       PRx: Calculated as the Pearson correlation coefficient between
%       MAP and ICP in a 300 s moving window. PRx is Fisher transformed
%       (artanh) to limit random correlation in the output. The first and
%       final 15 samples of PRx are set to NaN due to missing half of
%       window before and after start and end of signal, respectively.
% 
% Adjusted from Bart Formsma and Tahisa Robles
% Erasmus MC Sophia Childrens Hospital, Pediatric Intensive Care Unit in
% collaboration with Delft University of Technology, 3mE 
%
% Disclaimer: not validated for real time clinical use
% Eris van Twist, May 2024.

close all
clear all
clc

PRx = struct();
%% Calculations
% Fixed parameters: window and time sample (Ts = 1/Fs)
Ts = 10;                        % Seconds per sample
window = 300/(2*Ts);            % Sample points in 5 min window

PRxraw = [];
for nn = window+1:length(ICP)-window
    correlation = corrcoef(ICP(nn-window:nn+window), MAP(nn-window:nn+window)); 
    PRxraw = [PRxraw correlation(1,2)];
end

% Fisher transform
PRxfish = atanh(PRxraw);

% Finalize output with NaNs at start/end of PRx
PRx = [NaN(1,(150/Ts)), PRxfish, NaN(1,(150/Ts))];

%% Save to output
PRx.raw = PRxraw;
PRx.fisher_transform = PRxfish;
PRx.data = PRx

