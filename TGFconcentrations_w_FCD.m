function y0 = TGFconcentrations_w_FCD
%
%==========================================================================
% The best-fit parameters in Schmierer et al., 2008
%==========================================================================

p = csvread('parameters_beta.dat');
% y0 = zeros(22,1); %without mS2-YFP and without mS4 YFP [and with explicit
% receptors]
y0 = zeros(23,1);
k = p(1)/p(2);
y0(1) = p(9)*k/(1+k);        % cyt Smad2 
y0(2) = 0;                  %% cyt Smad2-Egfp
y0(3) = 0;           % cyt phospho Smad2
y0(4) = 0;           %% cyt phospho Smad2Egfp
y0(5) = p(10)/2;        % cyt Smad4
y0(6) = 0;           % cyt pSmad2-Smad4
y0(7) = 0;           %% cyt Smad2Egfp-Smad4
y0(8) = 0;           % cyt pSmad2-pSmad2
y0(9) = 0;           %% cyt Smad2Egfp-Smad2
y0(10) = 0;          %% cyt Smad2Egfp-Smad2Egfp
y0(11) = p(9)/(1+k);       % nuc Smad2
y0(12) = 0;       % nuc Smad2Egfp
y0(13) = 0;          % nuc phospho Smad2
y0(14) = 0;          %% nuc phospho Smad2Egfp
y0(15) = p(10)/2;    % nuc Smad4
y0(16) = 0;          % nuc Smad2-Smad4
y0(17) = 0;          %% nuc Smad2Egfp-Smad4
y0(18) = 0;          % nuc pSmad2-pSmad2
y0(19) = 0;          %% nuc Smad2Egfp-Smad2
y0(20) = 0;          %% nuc Smad2Egfp-Smad2Egfp
y0(21) = 1;%%%%          %% receptors *cytoplasmic
y0(22) = 0; %%         %% TGFbeta *cytoplasmic
y0(23) = 0;          %% active receptors
y0(24) = 1;        %IFFL intermediate
y0(25) = 1;        %IFFL output gene


% y0 = zeros(20,1); %with mS2-YFP and mS4 YFP
% k = p(1)/p(2);
% y0(1) = p(9)*k/(1+k);        % cyt Smad2 
% y0(2) = p(9)*k/(1+k);        % cyt Smad2-Egfp
% y0(3) = 0;           % cyt phospho Smad2
% y0(4) = 0;           % cyt phospho Smad2Egfp
% y0(5) = p(10)/2;        % cyt Smad4
% y0(6) = 0;           % cyt Smad2-Smad4
% y0(7) = 0;           % cyt Smad2Egfp-Smad4
% y0(8) = 0;           % cyt Smad2-Smad2
% y0(9) = 0;           % cyt Smad2Egfp-Smad2
% y0(10) = 0;          % cyt Smad2Egfp-Smad2Egfp
% y0(11) = p(9)/(1+k);       % nuc Smad2
% y0(12) = p(9)/(1+k);       % nuc Smad2Egfp
% y0(13) = 0;          % nuc phospho Smad2
% y0(14) = 0;          % nuc phospho Smad2Egfp
% y0(15) = p(10)/2;       % nuc Smad4
% y0(16) = 0;          % nuc Smad2-Smad4
% y0(17) = 0;          % nuc Smad2Egfp-Smad4
% y0(18) = 0;          % nuc Smad2-Smad2
% y0(19) = 0;          % nuc Smad2Egfp-Smad2
% y0(20) = 0;          % nuc Smad2Egfp-Smad2Egfp

% y0 = zeros(20,1); %without mS2-YFP and WITH mS4 YFP
% k = p(1)/p(2);
% y0(1) = p(9)*k/(1+k);        % cyt Smad2 
% y0(2) = 0;        % cyt Smad2-Egfp
% y0(3) = 0;           % cyt phospho Smad2
% y0(4) = 0;           % cyt phospho Smad2Egfp
% y0(5) = p(10)/2;        % cyt Smad4
% y0(6) = 0;           % cyt Smad2-Smad4
% y0(7) = 0;           % cyt Smad2Egfp-Smad4
% y0(8) = 0;           % cyt Smad2-Smad2
% y0(9) = 0;           % cyt Smad2Egfp-Smad2
% y0(10) = 0;          % cyt Smad2Egfp-Smad2Egfp
% y0(11) = p(9)/(1+k);       % nuc Smad2
% y0(12) = 0;       % nuc Smad2Egfp
% y0(13) = 0;          % nuc phospho Smad2
% y0(14) = 0;          % nuc phospho Smad2Egfp
% y0(15) = p(10)/2;       % nuc Smad4
% y0(16) = 0;          % nuc Smad2-Smad4
% y0(17) = 0;          % nuc Smad2Egfp-Smad4
% y0(18) = 0;          % nuc Smad2-Smad2
% y0(19) = 0;          % nuc Smad2Egfp-Smad2
% y0(20) = 0;          % nuc Smad2Egfp-Smad2Egfp

% y0 = zeros(20,1); %WITH mS2-YFP and without mS4 YFP
% k = p(1)/p(2);
% y0(1) = p(9)*k/(1+k);        % cyt Smad2 
% y0(2) = p(9)*k/(1+k);        % cyt Smad2-Egfp
% y0(3) = 0;           % cyt phospho Smad2
% y0(4) = 0;           % cyt phospho Smad2Egfp
% y0(5) = p(10)/2;        % cyt Smad4
% y0(6) = 0;           % cyt Smad2-Smad4
% y0(7) = 0;           % cyt Smad2Egfp-Smad4
% y0(8) = 0;           % cyt Smad2-Smad2
% y0(9) = 0;           % cyt Smad2Egfp-Smad2
% y0(10) = 0;          % cyt Smad2Egfp-Smad2Egfp
% y0(11) = p(9)/(1+k);       % nuc Smad2
% y0(12) = p(9)/(1+k);       % nuc Smad2Egfp
% y0(13) = 0;          % nuc phospho Smad2
% y0(14) = 0;          % nuc phospho Smad2Egfp
% y0(15) = p(10)/2;       % nuc Smad4
% y0(16) = 0;          % nuc Smad2-Smad4
% y0(17) = 0;          % nuc Smad2Egfp-Smad4
% y0(18) = 0;          % nuc Smad2-Smad2
% y0(19) = 0;          % nuc Smad2Egfp-Smad2
% y0(20) = 0;          % nuc Smad2Egfp-Smad2Egfp
% y0(21) = 1;
% y0(22) = 0;
% y0(23) = 0;



% y0(21) = 0.01;
% %Y
% y0(22) = 10;
%Z

% Y(:,11)+Y(:,13)+Y(:,16)+Y(:,18); nuclear Smad2
% Y(:,1)+Y(:,3)+Y(:,6)+Y(:,8); cytoplasmic Smad2
% Y(:,6)+Y(:,16); %total S24
% Y(:,5)+Y(:,6); %cytoplasmic S4
% Y(:,15)+Y(:,16); %nuclear S4




% Note that based on the koff and kon values, the lifetime of the Smad
% complex is ~1 minute. 