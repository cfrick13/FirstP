function p = TGFparameters_w_FCD
%

%==========================================================================
% The best-fit parameters in Schmierer et al., 2008
%==========================================================================


%%%%%%%%%%%%%%%%%%%%%
%test set of parameters
%%%%%%%%%%%%%%%%%%%%%%
p(1) = 0.0056;       % /s, kex, page 1 of Supp 0.0056;  
p(2) = 0.0026;       % /s, kin, page 1 of Supp 0.0026
p(3) = 0.000404;     % /nM-s, kphosp, Fig S4kex
p(4) = (0.016/8.7);    % /nM-s, kon = koff / Kdiss, Kdiss is in Fig S4
p(5) = 0.016;        % /s, koff, page 2 of Supp
p(6) = 5.7;          % CIF, Fig S4
p(7) = 0.00657;      % /nM-s, kdephosp, Fig S4
%increasing dephos rate increases robustness!!!
%and decreasesing decreases robustness
p(8) = 1;            % nM, PPase, Table S1 
p(9)= 89.1;             % nM, S2total           
p(10)= 101.6;        % nM, S4total       
p(11) = 0.074;          %/nM-s, kTGFbeta

p(12) = 10;     % kxy
p(13) = 10;     % kxz
p(14) = 10;     % kyz
p(15) = 0.0001;     % betaY
p(16) = 2;     % betaZ
p(17) = 0.0005;     % gammmaY;
p(18) = 100;     % gammaZ;
p(19) =1000; %differential binding


% 
% %%%%%%%%%%%%%%%%%%%%%
% %working FCD parameters
% %%%%%%%%%%%%%%%%%%%%%%
% p(1) = 0.0056;       % /s, kex, page 1 of Supp 0.0056;  
% p(2) = 0.0026;       % /s, kin, page 1 of Supp 0.0026
% p(3) = 0.000404;     % /nM-s, kphosp, Fig S4kex
% p(4) = (0.016/8.7);    % /nM-s, kon = koff / Kdiss, Kdiss is in Fig S4
% p(5) = 0.016;        % /s, koff, page 2 of Supp
% p(6) = 5.7;          % CIF, Fig S4
% p(7) = 0.00657;      % /nM-s, kdephosp, Fig S4
% %increasing dephos rate increases robustness!!!
% %and decreasesing decreases robustness
% p(8) = 1;            % nM, PPase, Table S1 
% p(9)= 89.1;             % nM, S2total           
% p(10)= 101.6;        % nM, S4total       
% p(11) = 0.074;          %/nM-s, kTGFbeta
% 
% p(12) = 10;     % kxy
% p(13) = 10;     % kxz
% p(14) = 10;     % kyz
% p(15) = 0.0001;     % betaY
% p(16) = 2;     % betaZ
% p(17) = 0.0005;     % gammmaY;
% p(18) = 100;     % gammaZ;
% p(19) =1000; %differential binding

