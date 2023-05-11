%% CONFIGURACION GLOBAL

%% OCTAVAS

octavas = [
       14.1    16      17.8
       17.8    20      22.4
       22.4    25      28.2
       28.2    31.5    35.5
       35.5    40      44.7
       44.7    50      56.2
       56.2    63      70.8
       70.8    80      89.1
       89.1   100     112  
      112     125     141  
      141     160     178  
      178     200     224  
      224     250     282  
      282     315     355  
      355     400     447  
      447     500     562  
      562     630     708  
      708     800     891  
      891    1000    1122  
     1122    1250    1413  
     1413    1600    1778  
     1778    2000    2239  
     2239    2500    2818  
     2818    3150    3548  
     3548    4000    4467  
     4467    5000    5623  
     5623    6300    7079  
     7079    8000    8913  
     8913   10000   11220  
%     11220   12220   14130  
%     14130   16000   17780  
%     17780   20000   22390  
];

% Potencia externa aplicada
Pext = zeros(size(octavas));

% Frecuencias centrales
c = octavas(:, 2);

Pext( (c >=  16 ) & (c <=  1000), 1 ) = 10   ;
Pext( (c >=  16 ) & (c <=  1000), 2 ) =  4.35;
Pext( (c >=  16 ) & (c <=  1000), 3 ) = 10   ;

Pext( (c >= 5000) & (c <= 10000), 1 ) = 100   ;
Pext( (c >= 5000) & (c <= 10000), 2 ) =  43.47;
Pext( (c >= 5000) & (c <= 10000), 3 ) = 100   ;

Pext( c == 1250, : ) = [ 20  8.70  20];
Pext( c == 1600, : ) = [ 35 15.20  35];
Pext( c == 2000, : ) = [ 50 21.74  50];
Pext( c == 2500, : ) = [ 80 36.96  80];
Pext( c == 3150, : ) = [100 39.13 100];
Pext( c == 4000, : ) = [150 45.65 150];

clear c



%% CONFIGURACION DE LOS MATERIALES

% Material aluminio
aluminio = struct();

% Densidad [kg/m^3], modulo de Young [Pa] y coeficiente de Poisson
aluminio.rho = 2700;
aluminio.E   = 70e9;
aluminio.v   = 0.33;

% Material aire
aire = struct();

% Densidad [kg/m^3], modulo de Young [Pa] y coeficiente de Poisson
aire.rho = 1.23;
aire.E   =    1;
aire.v   =    0;

% Velocidad del sonido [m/s]
aire.c = 343;



%% CONFIGURACION DE LAS SECCIONES

% Seccion de panel
panel = section(aluminio, 1, 1.25, 5e-3);
panel.k = 0.015;

% Seccion de aire
capa = section(aire, 1, 1.25, 50e-3);
capa.k = 0.010;


%% CONFIGURACION DE ESPACIOS

% Espacio de frecuencia
Nf = 10000-16;
f = linspace(16, 10000, Nf);

fprintf("df is %.3f\n", f(2) - f(1));


%% CALCULOS INICIALES

% Calcular las frecuencias fc y f11
panel.fc = fc(panel, capa);
panel.f11 = f11(panel, capa);


%% APARTADO 1: Numero de modos por banda

% Densidad n de los paneles y el aire
npanel = np(panel, (2 * pi) .* octavas(:, 2));
naire  = na(capa , (2 * pi) .* octavas(:, 2));

% Delta de frecuencia
dw = (2 * pi) .* (2^(1/6) - 2^(-1/6)) .* octavas(:, 2);

% Numero de modos
Npanel = npanel .* dw;
Naire  = naire  .* dw;

% Encontrar el punto donde ambos son al menos 5
[~, ip] = min( abs( Npanel - 5 ) );
[~, ia] = min( abs( Naire  - 5 ) );

i5 = max(ip, ia);
f5 = octavas(i5, 2);

figure;
loglog(octavas(:, 2), [Npanel, Naire], "-", 'LineWidth', 1.5);

xline(f5);

text(f5, 0.1, 0.1, sprintf("%.0f Hz", f5));

xlabel("Frecuencia [Hz]");
ylabel("Numero de nodos");
legend(["n_{p}", "n_{a}"]);


% SystemResponse(f, panel, capa, octavas, Pext);
% SystemResponse(linspace(16, 200, 10000), panel, capa, octavas, Pext);
% SystemResponse(linspace(2200, 2700, 10000), panel, capa, octavas, Pext);


%% APARTADO 2: Factor de perdidas de acoplamiento segun la frecuencia


% Array con los factores de amortiguamiento
eta = zeros(Nf, 2);

for i=1:Nf
    % Amortiguamiento de radiacion de potencia panel -> aire
    eta(i, 1) = kpa(panel, capa, 2 * pi * f(i));

    % Amortiguamiento de radiación de potencia aire -> panel
    eta(i, 2) = eta(i, 1) * np(panel, f(i)) / na(capa, 2 * pi * f(i));
end

DisplayK(f, eta(:, 1), eta(:, 2));

%% APARTADO 3: Ecuaciones de equilibrio de potencias

% Apartado teórico

%% APARTADO 4: Energia de los subsistemas

% Frecuencia y amortiguamiento bajo condicion de alta frecuencia
fhigh   = f( f >= f5 );
etahigh = eta( f >= f5, : );

% Energia del sistema
E = SystemEnergy(fhigh, panel, capa, etahigh, octavas, Pext);

% Dibujar la grafica
DisplayEnergy(fhigh, E);

% Energia del sistema
octhigh = octavas( octavas(:, 1) >= f5 , :);
etaoct = zeros(length(octhigh), 2);

for i=1:length(etaoct)
    [~, idx] = min( abs( f - octhigh(i, 2) ) );
    element = eta( idx, : );
    etaoct(i, :) = element;
end

Eoct = SystemEnergy(octhigh(:, 2), panel, capa, etaoct, octavas, Pext);

% Dibujar la grafica
DisplayEnergy(octhigh(:, 2), Eoct);


%% APARTADO 5: Velocidad media de los paneles y presion media del aire

% Para todo el espacio de frecuencias
% Velocidad media de los paneles
Vrms = sqrt( E(:, 1:2:5) ./ panel.m );

% Presión media de los paneles
Prms = sqrt( E(:, 2:2:4) .* (aire.rho * aire.c) ./ capa.V );

% Dibujar la grafica
DisplayRMS(fhigh, Vrms, Prms, octavas);

% Para las frecuencias de octavas
% Velocidad media de los paneles
Vrms = sqrt( Eoct(:, 1:2:5) ./ panel.m );

% Presión media de los paneles
Prms = sqrt( Eoct(:, 2:2:4) .* (aire.rho * aire.c) ./ capa.V );

% Dibujar la grafica
DisplayRMS(octhigh(:, 2), Vrms, Prms, octavas);

%% FUNCIONES DENSIDAD NODAL

% Densidad modal de una seccion panel
function N=np(panel, ~)
    N = sqrt(panel.material.rho * panel.Lz / panel.D) * panel.A / (4 * pi);
end

% Densidad modal de una seccion aire
function N=na(aire, w)
    % Dividir en partes para mayor legibilidad
    N1 = (aire.V / (pi * aire.material.c)) .* ((w ./ aire.material.c).^2);

    N2 = (w ./ aire.material.c) .* (aire.A / (4 * aire.material.c));

    N3 = aire.P / (8 * aire.material.c);

    N = N1 + N2 + N3;
end

% Densidad modal cruzada
function N=kpa(panel, aire, w)
    s = sigma(panel, aire, w);
    N = (panel.A * aire.material.rho * aire.material.c * s) / (panel.m * w);
end

%% FUNCIONES DE FRECUENCIAS DE RADIACION

% Calcula la frecuencua critica de la placa
function F = fc(panel, aire)
    F = ((aire.material.c ^ 2) / (2 * pi)) * sqrt((panel.m / panel.A) / panel.D);
end

% Calcula el primer modo superficial de la placa
function F = f11(panel, aire)
    F = ((aire.material.c ^ 2) / (4 * panel.fc)) * ((1 / (panel.Lx ^ 2)) + (1 / (panel.Ly ^ 2)));
end

% Eficiencia de radiacion
function S = sigma(panel, aire, w)
    % Frecuencia de excitacion
    f = w / (2 * pi);

    % Factor lambda
    l = sqrt( f / panel.fc );

    % Formula si f11 < fc / 2
    if panel.f11 < (panel.fc / 2)
        if f < panel.f11
            S = 4 * ((panel.P * (f / aire.material.c)) ^ 2);
        elseif f < panel.fc
            % Primer factor comun
            A = (aire.material.c * panel.P) / (panel.A * panel.fc);
            B = ((1 - (l^2)) * log((1+l) / (1-l))) + (2 * l);
            C = 4 * (pi ^ 2) * ((1 - (l ^ 2)) ^ (3 / 2));

            % Factor no comun
            D = 0;

            % Cambiar factor D si f < (fc / 2)
            if f < (panel.fc / 2)
                X = 2 * (( (2 * aire.material.c) / (panel.fc * (pi ^ 2)) ) ^ 2);
                Y = (1 - (2 * (l ^ 2))) / (panel.A * l * sqrt(1 - (l ^ 2)));

                D = X * Y;
            end

            S = (A * B / C) + D;
        else
            S = 1 / sqrt( 1 - (panel.fc / f) );
        end
    else
        if f < panel.fc
            S = 4 * ((panel.A * f / aire.material.c) ^ 2);
        else
            S = 1 / sqrt(1 - (f / panel.fc));
        end
    end
end



%% FUNCIONES DE CREACION DE SECCIONES

function S = section(material, Lx, Ly, Lz)
    % Inicializar struct
    S = struct();

    % Material
    S.material = material;

    % Longitudes de la seccion [m]
    S.Lx = Lx;
    S.Ly = Ly;
    S.Lz = Lz;
    
    % Area de la seccion [m^2]
    S.A = S.Lx * S.Ly;
    
    % Volumen de la seccion [m^3]
    S.V = S.A * S.Lz;

    % Perimetro de la seccion
    S.P = (2 * S.Lx) + (2 * S.Ly);
    
    % Masa de la seccion
    S.m = material.rho * S.V;

    % Rigidez de la seccion
    S.D = (material.E * (S.Lz ^ 3)) / (12 * (1 - (material.v ^ 2)));
end

%% FUNCIONES DE RESOLUCION

function E = SystemEnergy(f, panel, capa, eta, octavas, Pext)
    % Vector de energía
    E = zeros(length(f), 5);

    % Precalcular w
    w = (2 * pi) .* f;

    % Valores de la potencia
    fw = [1000 1250 1600 2000 2500 3150 4000 10000];
    fp = [
         10  4.35
         20  8.70
         35 15.20
         50 21.74
         80 36.96
        100 39.13
        150 45.65
        100 43.37
    ];

    for i=1:length(f)
        % Valores de la potencia
        P = fp( find(fw >= f(i), 1), :);

        % Amortiguamiento cruzado para el panel y el aire
        epa = eta(i, 1);
        eap = eta(i, 2);

        % Calcular la diagional
        D = diag([(panel.k + epa) (capa.k + (2 * eap)) (panel.k + (2 * epa)) (capa.k + (2 * eap)) (panel.k + epa)]);

        % Calcular el resto de la matriz de transferencia de potencia
        T = [
               0  -eap   0    0    0
             -epa   0  -epa   0    0
               0  -eap   0  -eap   0
               0    0  -epa   0  -epa
               0    0    0  -eap   0
        ];

        % Matrix de transferencia 
        M = D + T;
    
        % Vector de potencia introducida
        P = [P(1) 0 P(2) 0 P(1)] ./ w(i);

        % Resolver el sistema
        E(i, :) = M \ P';
    end
end

function SystemResponse(f, panel, capa, eta, octavas, Pext)
    % Vector de energia
    Erms = zeros(length(f), 5);

    % Precalcular w
    w = (2 * pi) .* f;
    
    for i=1:length(f)
        % Frecuencia central
        fhz = f(i);
    
        % Amortiguamiento cruzado para el panel y el aire
        kcp = eta(i, 1);
        kca = eta(i, 2);
    
        % Calcular la matriz de transferencia de potencia
        K = w(i) * [
            panel.k   -kcp      0       0       0
             -kca    capa.k   -kca      0       0
               0      -kcp   panel.k  -kcp      0
               0        0     -kca   capa.k   -kca
               0        0       0      -kcp  panel.k
        ];
    
        % Indice de octava
        io = 1;
    
        for j=1:length(octavas)
            if (fhz >= octavas(j, 1)) && (fhz < octavas(j, 3))
                io = j;
                break;
            end
        end
    
        % Vector de potencia introducida
        P = [Pext(io, 1) 0 Pext(io, 2) 0 Pext(io, 3)];
    
        % Resolver el sistema
        Erms(i, :) = K \ (P');
    end
    
    % Sacar la Vrms
    Vrms = log10( real( sqrt(Erms(:, 1:2:5) ./ panel.m) ) );
    
    % Sacar la Prms
    Prms = log10( real( sqrt(Erms(:, 2:2:4) .* (capa.material.rho * capa.material.c) ./ capa.V) ) );
    
    
    figure;
    hold on;
    
    xlabel("Frecuencia [Hz]");

    % Colores de los patches
    colores = ['b' 'y'];

    % Calcular maximo de Vrms.
    vmx = max(Vrms(Vrms ~= Inf & Vrms ~= -Inf));
    vmn = min(Vrms(Vrms ~= Inf & Vrms ~= -Inf));

    % Patch las octavas
    for j=1:length(octavas)
        % Comprobar si hay que dibujarlo
        if ~any( f >= octavas(j, 1) & f <= octavas(j, 2) )
            continue;
        end

        % Coger el minimo y maximo
        mn = octavas(j, 1);
        mx = octavas(j, 3);

        % Dibujar el patch
        patch([mn mn mx mx], [vmn vmx vmx vmn], 'k', 'FaceColor', colores(rem(j, 2) + 1), 'FaceAlpha', 0.1);

        % Dibujar la linea vertical del centro de la octava
    end

    yyaxis left;
    P1 = plot(f, Vrms(:, 1), '-' , 'Color', '#0072BD', 'LineWidth', 1.5);
    P2 = plot(f, Vrms(:, 2), '--', 'Color', '#0072BD', 'LineWidth', 1.5);
    P3 = plot(f, Vrms(:, 3), ':' , 'Color', '#0072BD', 'LineWidth', 2.5);
    
    ylabel("log_{10}(V_{RMS}) [m/s]");
    
    yyaxis right;
    A1 = plot(f, Prms(:, 1), '-' , 'Color', '#D95319', 'LineWidth', 1.5);
    A2 = plot(f, Prms(:, 2), '--', 'Color', '#D95319', 'LineWidth', 1.5);
    
    ylabel("log_{10}(P_{RMS}) [Pa]");

    legend([P1, P2, P3, A1, A2], ["Placa 1", "Placa 2", "Placa 3", "Capa 1", "Capa 2"]);
    
    hold off;
end

%% FUNCIONES GRAFICAS

function DisplayRMS(f, V, P, octavas)
    % Crear una nueva figura
    figure;
    hold on;

    % Colores de los patches
    colores = ['b' 'y'];

    % Calcular maximo de Vrms.
    vmx = max(V(V ~= Inf & V ~= -Inf));
    vmn = min(V(V ~= Inf & V ~= -Inf));

    % Patch las octavas
    for j=1:length(octavas)
        % Comprobar si hay que dibujarlo
        if ~any( f >= octavas(j, 1) & f <= octavas(j, 2) )
            continue;
        end

        % Coger el minimo y maximo
        mn = octavas(j, 1);
        mx = octavas(j, 3);

        % Dibujar el patch
        patch([mn mn mx mx], [vmn vmx vmx vmn], 'k', 'FaceColor', colores(rem(j, 2) + 1), 'FaceAlpha', 0.1);
    end

    % Dibujar la Vrms
    yyaxis left;
    P1 = plot(f, V(:, 1), '-' , 'Color', '#0072BD', 'LineWidth', 1.5);
    P2 = plot(f, V(:, 2), '--', 'Color', '#0072BD', 'LineWidth', 1.5);
    P3 = plot(f, V(:, 3), ':' , 'Color', '#0072BD', 'LineWidth', 2.5);
    
    ylabel("V_{RMS} [m/s]");

    % Dibujar la Prms
    yyaxis right;
    A1 = plot(f, P(:, 1), '-' , 'Color', '#D95319', 'LineWidth', 1.5);
    A2 = plot(f, P(:, 2), '--', 'Color', '#D95319', 'LineWidth', 1.5);
    
    ylabel("P_{RMS} [Pa]");

    % Leyenda y Xlabel
    xlabel("Frecuencia [Hz]");
    legend([P1, P2, P3, A1, A2], ["Placa 1", "Placa 2", "Placa 3", "Capa de aire 1", "Capa de aire 2"]);
    
    hold off;
end

function DisplayEnergy(f, E)
    % Crear una nueva figura
    figure;
    hold on;

    % Dibujar las graficas
    plot(f, E(:, 1), ':' , 'LineWidth', 1.5);
    plot(f, E(:, 3), '--', 'LineWidth', 1.5);
    plot(f, E(:, 5), '-' , 'LineWidth', 1.5);
    plot(f, E(:, 2), ':' , 'LineWidth', 1.5);
    plot(f, E(:, 4), '-' , 'LineWidth', 1.5);

    % Y axis as log
    %set(gca, 'YScale', 'log');

    % Labels y leyenda
    xlabel("Frecuencia [Hz]");
    ylabel("Energia [J]");

    legend(["Placa 1", "Placa 2", "Placa 3", "Capa de aire 1", "Capa de aire 2"]);

    hold off
end

function DisplayK(f, kpa, kap)
    % Crear una nueva figura
    figure;
    hold on;

    plot(f, log10(kpa),  '-', 'LineWidth', 1.5);
    plot(f, log10(kap), '--', 'LineWidth', 1.5);

    % Labels
    xlabel("Frecuencia [Hz]");
    ylabel("log_{10}(\eta)");

    % Leyenda
    legend(["\eta_{pa}", "\eta_{ap}"]);

    hold off;
end
