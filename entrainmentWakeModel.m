function [Uw, Dw, Uw_lin, Dw_lin, kE] = entrainmentWakeModel(x, Ct, E)
%
% [Uw, Dw] = entrainmentWakeModel(x, Ct, E)
%
% where:
% x is the distance downstream of the rotor, normalized by rotor diameter;
% Ct is the thrust coefficient;
% E is the entrainment coefficient (if not supplied, E = 0.15 is assumed);
% Uw is the velocity in the wake, normalized by upstream velocity;
% Dw is the diameter of the wake, normalized by rotor diameter.
%
% To produce sample output and plots, run without any input.
%
% References:
%
% Morton, B. R. 1961. ?On a Momentum-Mass Flux Diagram for Turbulent Jets, 
% Plumes and Wakes.? Journal of Fluid Mechanics 10 (1): 101?12.
%
% Luzzatto-Fegiz, P. 2018 ""A one-parameter model for turbine wakes from 
% the entrainment hypothesis" Journal of Physics: Conf. Series 1037:072019
%
% Paolo Luzzatto-Fegiz, May 2018

if ~exist('Ct', 'var')
    Ct = .8;
end
if ~exist('E', 'var')
    E = 0.15;
end
if ~exist('x', 'var')
    x = linspace(0,20,100);
end

X = 6*E*(2/Ct)^0.5 * x + (1-Ct)^(3/4) / (1-(1-Ct)^(1/2))^(3/2);
% x_v =  - (1-Ct)^(3/4) / (1-(1-Ct)^(1/2))^(3/2) /( 6*E*(2/Ct)^0.5 )

Uw = X.^(2/3) ./ (X.^(2/3) + 1) ;
Dw = sqrt(Ct/2) * (X.^(2/3) + 1 )./ X.^(1/3);

if nargout > 2
    alpha = 6*E*(2/Ct)^0.5;
    beta = (1-Ct)^(3/4) / (1-(1-Ct)^(1/2))^(3/2);
    [kE, x_r] = E_Ct2k(E, Ct);
    betaTilde = alpha*x_r + beta;
    DiE = sqrt(Ct/2)*( (betaTilde^(2/3)+1)/betaTilde^(1/3) ...
        - alpha*x_r*(betaTilde^(2/3)-1)/3/betaTilde^(4/3) );
    Dw_lin = DiE + 2*kE*x;
    Uw_lin = 0.5*(1-sqrt(1-2*Ct./Dw_lin.^2 ) );
end

if nargout == 0
    figure
    subplot(2,1,1)
    plot(x, Uw,'r');hold on
    title('Entrainment model')
    ylabel('U_w/U_o')
    subplot(2,1,2)
    plot(x, Dw, 'r')
    hold on
    ylabel('D_w/D')
    xlabel('x/D')
    
    
    % small-x approximation
    alpha = 6*E*(2/Ct)^0.5;
    beta = (1-Ct)^(3/4) / (1-(1-Ct)^(1/2))^(3/2);
    
    % linear expansion 
    [kE, x_r] = E_Ct2k(E, Ct);
    betaTilde = alpha*x_r + beta;
    DiE = sqrt(Ct/2)*( (betaTilde^(2/3)+1)/betaTilde^(1/3) ...
        - alpha*x_r*(betaTilde^(2/3)-1)/3/betaTilde^(4/3) );
    DwE_lin = DiE + 2*kE*x;
%     plot(x, DwE_lin, 'r--')
    axis([0 10 1 2])
    
    subplot(2,1,1)
    a=0.5*(1-sqrt(1-Ct));
    UwE_lin = 0.5*(1+sqrt(1-2*Ct./DwE_lin.^2 ) );
%     plot(x, UwE_lin,'r--');hold on
    
    
%     % compare momentum conservation to park model
%     figure
%     k = 0.06; % spreading paramer
%     a = 0.5*(1-sqrt(1-Ct));
%     Di = sqrt((1 - a)/(1 - 2*a));
%     Dpark = Di + 2*k*x;
%     Upark = 1 - 2*a./(1+2*k*x./Di).^2;
%     thrust = .5*Ct*pi/4; % nondimensionalized by upstream vel., rotor diam.
%     momentumDeficitEntrainment = pi/4*Dw.^2 .* Uw.*(1-Uw);
%     momentumDeficitPark = pi/4*Dpark.^2 .* Upark.*(1-Upark);
%     plot(x,x*0+thrust)
%     hold on
%     plot(x,momentumDeficitEntrainment,'r.:')
%     plot(x,momentumDeficitPark)
%     legend('Thrust','Wake momentum deficit, entrainment model',...
%         'Wake momentum deficit, park model','location','nw')
%     xlabel('x/D')
    

end
