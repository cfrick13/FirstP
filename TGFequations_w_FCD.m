function lhs = TGFequations_w_FCD(t,y,dy)
%constitutively nuclear Smad4 is modeled with a really high import rate for
%cytoplasmic Smad4 to mimic the nuclear import in cell

%
%==========================================================================
% Evaluate the LHS of the dimensional ODE's describing TGF-beta signalling.
% Based on Schmierer et al., 2008
%==========================================================================

global Tgfon Tgfoff Tgfbasal Tgfz Time  Basal

% Evaluating the left hand sides of the ODE's written in standard implicit
% for, f(t,y,y')=0

p = load('parameters_beta.dat');
lhs = zeros(23,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STANDARD MODEL (from paper)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lhs(1)  = -dy(1)  + p(1)*y(11) - p(2)*y(1) - p(3)*y(23)*y(1);  
lhs(2)  = -dy(2)  + p(1)*y(12) - p(2)*y(2) - p(3)*y(23)*y(2);
lhs(3)  = -dy(3)  + p(1)*y(13) - p(2)*y(3) + p(3)*y(23)*y(1) - p(4)*y(3)*(y(5) + 2*y(3) + y(4)) + p(5)*(y(6) + 2*y(8) + y(9));
lhs(4)  = -dy(4)  + p(1)*y(14) - p(2)*y(4) + p(3)*y(23)*y(2) - p(4)*y(4)*(y(5) + y(3) + 2*y(4)) + p(5)*(y(7) + y(9) + 2*y(10));
lhs(5)  = -dy(5)  + p(2)*y(15) - p(2)*y(5) - p(4)*y(5)*(y(3) + y(4)) + p(5)*(y(6)+y(7));
lhs(6)  = -dy(6)  + p(4)*y(3)*y(5) - p(5)*y(6) - p(2)*p(6)*y(6);
lhs(7)  = -dy(7)  + p(4)*y(4)*y(5) - p(5)*y(7) - p(2)*p(6)*y(7);
lhs(8)  = -dy(8)  + p(4)*y(3)^2 - p(5)*y(8) - p(2)*p(6)*y(8);
lhs(9)  = -dy(9)  + p(4)*y(4)*y(3) - p(5)*y(9) - p(2)*p(6)*y(9);
lhs(10) = -dy(10) + p(4)*y(4)^2 - p(5)*y(10) - p(2)*p(6)*y(10);
lhs(11) = -dy(11) + p(2)*y(1) - p(1)*y(11) + p(7)*p(8)*y(13);
lhs(12) = -dy(12) + p(2)*y(2) - p(1)*y(12) + p(7)*p(8)*y(14);
lhs(13) = -dy(13) + p(2)*y(3) - p(1)*y(13) - p(7)*p(8)*y(13) - p(4)*y(13)*(y(15) + 2*y(13) + y(14)) + p(5)*(y(16) + 2*y(18) + y(19));
lhs(14) = -dy(14) + p(2)*y(4) - p(1)*y(14) - p(7)*p(8)*y(14) - p(4)*y(14)*(y(15) + y(13) + 2*y(14)) + p(5)*(y(17) + y(19) + 2*y(20));
lhs(15) = -dy(15) + p(2)*y(5) - p(2)*y(15) - p(4)*y(15)*(y(13) + y(14)) + p(5)*(y(16) + y(17));
lhs(16) = -dy(16) + p(4)*y(13)*y(15) - p(5)*y(16) + p(2)*p(6)*y(6);
lhs(17) = -dy(17) + p(4)*y(14)*y(15) - p(5)*y(17) + p(2)*p(6)*y(7);
lhs(18) = -dy(18) + p(4)*y(13)^2 - p(5)*y(18) + p(2)*p(6)*y(8);
lhs(19) = -dy(19) + p(4)*y(14)*y(13) - p(5)*y(19) + p(2)*p(6)*y(9);
lhs(20) = -dy(20) + p(4)*y(14)^2 - p(5)*y(20) + p(2)*p(6)*y(10);

lhs(21) = -dy(21) - p(11)*y(21)*y(22);  %Receptors
lhs(22) = -dy(22) - p(11)*y(21)*y(22); %TGF
lhs(23) = -dy(23) + p(11)*y(21)*y(22); %Active Receptors

% lhs(24) = -dy(24) + p(15)*(y(16)/p(12))/(1 + y(16)/p(12)) - p(17)*y(24);  %IFFL intermediate
% % lhs(24) = -dy(24) + p(15)*(y(16)/p(12))/(1 + y(16)/p(12)) + p(17)*y(24);  %IFFL intermediate
% lhs(25) = -dy(25) + p(16)*(y(16)/p(13))/((1 + y(16)/p(13))*(1/(1 + y(24)/p(14)))) - p(18)*y(25);  %IFFL output gene


%working FCD circuit equations
lhs(24) = -dy(24) + p(15)*(y(16)) - p(17)*y(24);  %IFFL intermediate
% lhs(24) = -dy(24) + p(15)*(y(16)/p(12))/(1 + y(16)/p(12)) + p(17)*y(24);  %IFFL intermediate
lhs(25) = -dy(25) + p(16)*(y(16)/(y(24)/p(19))) - p(18)*y(25);


















