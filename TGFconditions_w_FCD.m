function c = TGFconditions_w_FCD

%
%==========================================================================
% The best-fit parameters in Schmierer et al., 2008
%==========================================================================

%% standard
c = zeros(13,1);

c(1)    =   8000*60  ;   %tn
% c(2)    =   0.1   ;   %Tgfoff
% c(3)    =   0.1   ;   %Tgfbasal
% c(4)    =   1    ;  %Tgfon
c(2)    =   0.02   ;   %Tgfoff
c(3)    =   0.02 ;   %Tgfbasal
c(4)    =   0.2  ;  %Tgfon
c(5)    =   c(1) ;    %time
c(6)    =   7000*60   ;   %pulse
c(7)    =   100      ;    %omega, (2pi/omega*60)
c(8)    =   0*60     ;    %Decay
c(9)    =   2       ;    %Basal time multiplier
c(10)   =   0.5     ;    %DC (fraction of TGF)
c(11)   =   0.01     ;   %Amplitude (fraction of TGF)
c(12)   =   400*60  ;    %RampPulse
c(13)   =   1      ;    %Parameter Variation

% c(2)    =   0.019      %Tgfoff
% c(3)    =   0.019      %Tgfbasal

%% non standard
% c = zeros(13,1);
% 
% c(1)    =   8000*60  ;   %tn
% c(2)    =   0.037   ;   %Tgfoff
% c(3)    =   0.037   ;   %Tgfbasal
% c(4)    =   0.2     ;  %Tgfon
% c(5)    =   1000*60  ;    %time
% c(6)    =   7000*60   ;   %pulse
% c(7)    =   100      ;    %omega, (2pi/omega*60)
% c(8)    =   0*60     ;    %Decay
% c(9)    =   2       ;    %Basal
% c(10)   =   0.5     ;    %DC (fraction of TGF)
% c(11)   =   0.01     ;   %Amplitude (fraction of TGF)
% c(12)   =   400*60  ;    %RampPulse
% c(13)   =   1      ;    %Parameter Variation
