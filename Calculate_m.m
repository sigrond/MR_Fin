function [ m ] = Calculate_m(Temperature,lambda, liquid)

% lambda in nm
% Temperature in deg Celsius

T_ref = 20;

switch liquid

case 'EG'
A0 = 1.42522;
% Cauchy dispersion parameters:
A1 = 1984.54821;
A2 = 104174632;
m = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case '2EG'
A0 = 1.4406; 
% Cauchy dispersion parameters:
A1 = 1984.54821;
A2 = 104174632;
m = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case '3EG'
A0 = 1.4492;
% Cauchy dispersion parameters:
A1 = 1984.54821;
A2 = 104174632;
m = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case '4EG'
A0 = 1.4511;
% Cauchy dispersion parameters:
A1 = 1984.54821;
A2 = 104174632;
m = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case 'DMI'
A0 = 1.4664;  
% Cauchy dispersion parameters:
A1 = 1984.54821;
A2 = 104174632;
m = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case 'DMSO'
A0 = 1.4585;  
% Cauchy dispersion parameters:
A1 = 8584;
A2 = -485600000;
m = -0.00044*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case 'Glycerine'
A0 = 1.45797;  
% Cauchy dispersion parameters:
A1 = 5980;
A2 = -360000000;
    
m = -0.00044*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;
    
case 'H2O'
% dispersion parameters for water:
a0 = 0.244257733; 
a1 = 9.74634476*10^(-3);
a2 = (-3.73234996)*10^(-3); 
a3 = 2.68678472*10^(-4);
a4 = 1.5892057*10^(-3);
a5 = 2.45934259*10^(-3);
a6 = 0.90070492;
a7 = (-1.66626219)*10^(-2);
lambda_uv = 0.229202;
lambda_IR = 5.432937;
laser = lambda/589;

T_K = 273.15 + Temperature;

if ((0<=Temperature)&&(Temperature<=40))
    % Thiesen formula for 0 - 40 degC from Tanaka
    b1 = -3.983035;
    b2 = 301.797;
    b3 = 522528.9;
    b4 = 69.34881;
%     b5 = 999.974950; % Standard Mean Ocean Water (isotopic composition)
    b5 = 999.972; % tap water (Chappuis)
    Rho = ...
        b5 * ( 1 - ((Temperature + b1)^2)*(Temperature + b2)/(b3 * (Temperature + b4)));
    Rho = Rho/1000;
else
    % from McCutcheon book
    Rho = 1 - (Temperature + 288.9414)/(508929.2*((Temperature + 68.12963))) ...
        *((Temperature - 3.9863)^2);
%     Rho*1000
end    
F = Rho * (a0 + a1*Rho + a2*T_K/273.15 + a3*laser^2*T_K/273.15 + a4/laser^2 ...
    + a5/(laser^2 - lambda_uv^2) + a6/(laser^2 - lambda_IR^2) + a7*Rho^2);
m = sqrt((2*F+1)/(1-F));

% added by Tho Do Duc
    case 'Meth'
        % Cauchy dispersion parameters:
        A0=1.3195;
        A1=3053.64419;
        A2=-34163639.3011;
        A3=2.622128e+12;
        m=-0.00040*(Temperature-T_ref)+A0+A1/(lambda^2)+A2/(lambda^4)+A3/(lambda^6);
        

    case 'EG+6G'
        A0 = 1.4257;
        % Cauchy dispersion parameters EG+temperaturowa zaleznosc roztworu EG+6G:
        A1 = 1984.54821;
        A2 = 104174632; 
        m = -0.0003*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;
    case 'H2O+2EG'
        % dispersion parameters for water:
        a0 = 0.244257733;
        a1 = 9.74634476*10^(-3);
        a2 = (-3.73234996)*10^(-3);
        a3 = 2.68678472*10^(-4);
        a4 = 1.5892057*10^(-3);
        a5 = 2.45934259*10^(-3);
        a6 = 0.90070492;
        a7 = (-1.66626219)*10^(-2);
        lambda_uv = 0.229202;
        lambda_IR = 5.432937;
        laser = lambda/589;

        T_K = 273.15 + Temperature;

        if ((0<=Temperature)&&(Temperature<=40))
            % Thiesen formula for 0 - 40 degC from Tanaka
            b1 = -3.983035;
            b2 = 301.797;
            b3 = 522528.9;
            b4 = 69.34881;
            %     b5 = 999.974950; % Standard Mean Ocean Water (isotopic composition)
            b5 = 999.972; % tap water (Chappuis)
            Rho = ...
                b5 * ( 1 - ((Temperature + b1)^2)*(Temperature + b2)/(b3 * (Temperature + b4)));
            Rho = Rho/1000;
        else
            % from McCutcheon book
            Rho = 1 - (Temperature + 288.9414)/(508929.2*((Temperature + 68.12963))) ...
                *((Temperature - 3.9863)^2);
            %     Rho*1000
        end
        F = Rho * (a0 + a1*Rho + a2*T_K/273.15 + a3*laser^2*T_K/273.15 + a4/laser^2 ...
            + a5/(laser^2 - lambda_uv^2) + a6/(laser^2 - lambda_IR^2) + a7*Rho^2);
        mh = sqrt((2*F+1)/(1-F));
        A0 = 1.4406;
        % Cauchy dispersion parameters:
        A1 = 1984.54821;
        A2 = 104174632;
        mEG = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;
        m = (mh + mEG)/2;
       case 'H2O+3EG'
           A0 = 1.4492;
           % Cauchy dispersion parameters:
           A1 = 1984.54821;
           A2 = 104174632;
           m3EG = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;
           
           % dispersion parameters for water:
        a0 = 0.244257733;
        a1 = 9.74634476*10^(-3);
        a2 = (-3.73234996)*10^(-3);
        a3 = 2.68678472*10^(-4);
        a4 = 1.5892057*10^(-3);
        a5 = 2.45934259*10^(-3);
        a6 = 0.90070492;
        a7 = (-1.66626219)*10^(-2);
        lambda_uv = 0.229202;
        lambda_IR = 5.432937;
        laser = lambda/589;

        T_K = 273.15 + Temperature;

        if ((0<=Temperature)&&(Temperature<=40))
            % Thiesen formula for 0 - 40 degC from Tanaka
            b1 = -3.983035;
            b2 = 301.797;
            b3 = 522528.9;
            b4 = 69.34881;
            %     b5 = 999.974950; % Standard Mean Ocean Water (isotopic composition)
            b5 = 999.972; % tap water (Chappuis)
            Rho = ...
                b5 * ( 1 - ((Temperature + b1)^2)*(Temperature + b2)/(b3 * (Temperature + b4)));
            Rho = Rho/1000;
        else
            % from McCutcheon book
            Rho = 1 - (Temperature + 288.9414)/(508929.2*((Temperature + 68.12963))) ...
                *((Temperature - 3.9863)^2);
            %     Rho*1000
        end
        F = Rho * (a0 + a1*Rho + a2*T_K/273.15 + a3*laser^2*T_K/273.15 + a4/laser^2 ...
            + a5/(laser^2 - lambda_uv^2) + a6/(laser^2 - lambda_IR^2) + a7*Rho^2);
        mh = sqrt((2*F+1)/(1-F));
        m = ( mh + m3EG ) / 2;
    case 'Au'
        if lambda == 532.07
            m = complex( 0.46636,2.40877);
        elseif lambda == 654.25
            m = complex( 0.16573,3.17963);
        end
%         c = 299792458;
%         h = 6.6260755e-34;
%         e = 1.60217733e-19;
%         omega_pl = 9.096; % eV;
%         gamma = 0.072; % eV;
%         omega = c * h / ( e * lambda * 1e-9 );
%         eps_infty = 9.84;
%         eps = eps_infty - ( omega_pl ^ 2 ) / ( omega^2 + complex(0,1) * gamma * omega );
%         m = sqrt( eps );
        case 'Au r < 20[ nm ]'
        c = 299792458;
        h = 6.6260755e-34;
        e = 1.60217733e-19;
        omega_pl = 9.096;% eV;
        gamma = 0.072;% eV;
        V_fermi = 1.4e6;%(M/S), E8(SM/S)
        r = 10; %[nm]                       TO  DO !!!!!!!!!!!!!!!
        gam = gamma + C*(V_fermi/(r*1e-9))*h/(2*pi*e);
        omega = c*h / ( e*lambda*1e-9 );
        eps_infty = 9.84;
        eps = eps_infty - ( omega_pl ^ 2 ) / ( omega^2 + complex(0,1) * gamma * omega );
    otherwise
        m = 0;
           s = sprintf('Liquid not found \nm = %1.f !!!',m)
       errordlg(s,'m Error');     
end




