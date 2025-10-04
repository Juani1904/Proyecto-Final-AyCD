% Level-2 MATLAB S-Function:
% S-function (System-function): La S significa que es una función de Simulink.
% Y el Nivel 2 simplemente indica que esta función se escribe con el lenguaje de
% programación: MATLAB que posee objetos con métodos (setup, Outputs, Update, etc.),
% en lugar de ser una Función de MATLAB simple (function y = ...).

function Simulacion(block)
    %MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
    %   The MATLAB S-function is written as a MATLAB function with the
    %   same name as the S-function. Replace 'msfuntmpl_basic' with the
    %   name of your S-function.
    %   Copyright 2003-2018 The MathWorks, Inc.

    %%
    %% The setup method is used to set up the basic attributes of the
    %% S-function such as ports, parameters, etc. Do not add any other
    %% calls to the main body of the function.
    %%
    setup(block);
end %function


%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C MEX counterpart: mdlInitializeSizes
function setup(block)
    
    %Primeramente determinamos el numero de puertos del bloque
    %INPUTS
    block.NumInputPorts  = 4;
    %OUTPUTS
    block.NumOutputPorts = 0;

    % Setup port properties to be inherited or dynamic:
    % - INHERITED: Las propiedades del puerto se ajustan estrictamente en función de las señales conectadas.
    %block.SetPreCompInpPortInfoToInherited;
    % - DYNAMIC: Las propiedades del puerto pueden ser establecidas dinámicamente en tiempo de ejecución. Como
    % por ejemplo: block.InputPort(1).Dimensions.
    block.SetPreCompInpPortInfoToDynamic;
    % block.SetPreCompOutPortInfoToDynamic;


    % Override input port properties:
    %En principio se trabajara con propiedades heredadas, no se va a
    %castear aca desde la funcion
    block.InputPort(1).Dimensions  = 1;  %[1x1]
    block.InputPort(1).DatatypeID  = 0;      % double
    block.InputPort(1).Complexity  = 'Real';
    block.InputPort(1).SamplingMode = 'Sample'; %(Propiedad agregada al código base)
    block.InputPort(1).DirectFeedthrough = true;

    block.InputPort(2).Dimensions  = 1; %[1x1]
    block.InputPort(2).DatatypeID  = 0;      % double
    block.InputPort(2).Complexity  = 'Real';
    block.InputPort(2).SamplingMode = 'Sample';
    block.InputPort(2).DirectFeedthrough = true;

    block.InputPort(3).Dimensions  = 1; %[1x1]
    block.InputPort(3).DatatypeID  = 0;      % double
    block.InputPort(3).Complexity  = 'Real';
    block.InputPort(3).SamplingMode = 'Sample';
    block.InputPort(3).DirectFeedthrough = true;

    block.InputPort(4).Dimensions  = 1; %[1x1]
    block.InputPort(4).DatatypeID  = 0;      % double
    block.InputPort(4).Complexity  = 'Real';
    block.InputPort(4).SamplingMode = 'Sample';
    block.InputPort(4).DirectFeedthrough = true;
    % 
    % block.InputPort(5).Dimensions  = 1; %[1x1]
    % block.InputPort(5).DatatypeID  = 0;      % double
    % block.InputPort(5).Complexity  = 'Real';
    % block.InputPort(5).SamplingMode = 'Sample';
    % block.InputPort(5).DirectFeedthrough = true;
    % 
    % block.InputPort(6).Dimensions  = 1; %[1x1] % Actuador_TWISTLOCKS
    % block.InputPort(6).DatatypeID  = 0;      % double
    % block.InputPort(6).Complexity  = 'Real';
    % block.InputPort(6).SamplingMode = 'Sample';
    % block.InputPort(6).DirectFeedthrough = true;
    % 
    % block.InputPort(7).Dimensions  = [1 24]; % N_containers
    % block.InputPort(7).DatatypeID  = 0;      % double
    % block.InputPort(7).Complexity  = 'Real';
    % block.InputPort(7).SamplingMode = 'Sample';
    % block.InputPort(7).DirectFeedthrough = true;
    % 
    % block.InputPort(8).Dimensions  = [1 22]; % AlturasLineasTransicion
    % block.InputPort(8).DatatypeID  = 0;      % double
    % block.InputPort(8).Complexity  = 'Real';
    % block.InputPort(8).SamplingMode = 'Sample';
    % block.InputPort(8).DirectFeedthrough = true;
    % 
    % block.InputPort(9).Dimensions  = 1; %[1x1] % MANUAL-AUTOMATICO
    % block.InputPort(9).DatatypeID  = 8;        % boolean
    % block.InputPort(9).Complexity  = 'Real';
    % block.InputPort(9).SamplingMode = 'Sample';
    % block.InputPort(9).DirectFeedthrough = true;
    % 
    % block.InputPort(10).Dimensions  = 1; %[1x1] % CONTROL BALANCEO
    % block.InputPort(10).DatatypeID  = 8;        % boolean
    % block.InputPort(10).Complexity  = 'Real';
    % block.InputPort(10).SamplingMode = 'Sample';
    % block.InputPort(10).DirectFeedthrough = true;
    % 
    % block.InputPort(11).Dimensions  = 1; %[1x1] % EMERGENCIA
    % block.InputPort(11).DatatypeID  = 0;      % double
    % block.InputPort(11).Complexity  = 'Real';
    % block.InputPort(11).SamplingMode = 'Sample';
    % block.InputPort(11).DirectFeedthrough = true;
    % 
    % block.InputPort(12).Dimensions  = 1; %[1x1] % PESO NOMINAL/SOBRECARGA
    % block.InputPort(12).DatatypeID  = 0;      % double
    % block.InputPort(12).Complexity  = 'Real';
    % block.InputPort(12).SamplingMode = 'Sample';
    % block.InputPort(12).DirectFeedthrough = true;
    % 
    % block.InputPort(13).Dimensions  = 1; %[1x1] % Masa estimada.
    % block.InputPort(13).DatatypeID  = 0;      % double
    % block.InputPort(13).Complexity  = 'Real';
    % block.InputPort(13).SamplingMode = 'Sample';
    % block.InputPort(13).DirectFeedthrough = true;

    % % Register parameters:
    %Los parametros se ingresan como valores fijos desde el bloque de la
    %Sfunction, se pueden usar para setear valores fijos. Por ahora lo
    %dejamos comentado
    block.NumDialogPrms     = 12;
    block.DialogPrmsTunable = repmat({'Nontunable'}, 1, block.NumDialogPrms);
    %{
    %   Dialog parameters can be either tunable or nontunable. A tunable parameter is a parameter that a
    %   user can change while the simulation is running. For example, if you have three dialog box
    %   parameters, you set them like the following line:
    %   block.DialogPrmsTunable = {'Tunable','Nontunable','Nontunable'};
    %   https://la.mathworks.com/help/simulink/sfg/dialog-parameters-matlab.html
    %}

    % Register sample times:
    %  [0 offset]            : Continuous sample time
    %  [positive_num offset] : Discrete sample time
    %  [-1, 0]               : Inherited sample time
    %  [-2, 0]               : Variable sample time
    block.SampleTimes = [0.1 0]; %(Sets the S-Function sample time). %0.01

    % Specify the block simStateCompliance. The allowed values are:
    %    'UnknownSimState', < The default setting; warn and assume DefaultSimState
    %    'DefaultSimState', < Same sim state as a built-in block
    %    'HasNoSimState',   < No sim state
    %    'CustomSimState',  < Has GetSimState and SetSimState methods
    %    'DisallowSimState' < Error out when saving or restoring the model sim state
    block.SimStateCompliance = 'DefaultSimState';

    %% -----------------------------------------------------------------
    %% The MATLAB S-function uses an internal registry for all
    %% block methods. You should register all relevant methods
    %% (optional and required) as illustrated below. You may choose
    %% any suitable name for the methods and implement these methods
    %% as local functions within the same file. See comments
    %% provided for each function for more information.
    %% -----------------------------------------------------------------

    block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
    % block.RegBlockMethod('InitializeConditions', @InitializeConditions);
    block.RegBlockMethod('Start', @Start);
    block.RegBlockMethod('Outputs', @Outputs);     % Required
    % block.RegBlockMethod('Update', @Update);
    % block.RegBlockMethod('Derivatives', @Derivatives);
    block.RegBlockMethod('Terminate', @Terminate); % Required
end %setup

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C MEX counterpart: mdlSetWorkWidths
function DoPostPropSetup(block)
    block.NumDworks = 1;

    block.Dwork(1).Name            = 'initialized';
    block.Dwork(1).Dimensions      = 1;
    block.Dwork(1).DatatypeID      = 8;      % boolean
    block.Dwork(1).Complexity      = 'Real'; % real
    block.Dwork(1).UsedAsDiscState = true;

%     block.Dwork(2).Name            = 'alto_headblock';
%     block.Dwork(2).Dimensions      = 1;
%     block.Dwork(2).DatatypeID      = 0;      % int32 (integer type)
%     block.Dwork(2).Complexity      = 'Real'; % real
%     block.Dwork(2).UsedAsDiscState = true;
end %DoPostPropSetup

%%
%% Start:
%%   Functionality    : Called once at start of model execution. If you
%%                      have states that should be initialized once, this
%%                      is the place to do it.
%%   Required         : No
%%   C MEX counterpart: mdlStart
function Start(block)
  close all;
  persistent fig img;

  initialized = false;

  x_t=block.InputPort(1).Data;
  x_l=block.InputPort(2).Data;
  y_l=block.InputPort(3).Data;
  %clock=block.InputPort(13).Data;
  %m_l_estimada = block.InputPort(5).Data;

  W_c = block.DialogPrm(1).Data; % Ancho de container (estándar): 2.438 m.
  H_c = block.DialogPrm(2).Data; % Alto de container (estándar): 2.591 m.
  alto_carro = block.DialogPrm(3).Data; % 1.500 m.
  alto_headblock = block.DialogPrm(4).Data; % 0.850 m.
  alto_spreader = block.DialogPrm(5).Data; % 0.750 m.
  ancho_carro = block.DialogPrm(6).Data; % 7.500 m.
  ancho_headblock = block.DialogPrm(7).Data; % 2.100 m.
  ancho_spreader = block.DialogPrm(8).Data; % 2.500 m.
  color_carro = block.DialogPrm(9).Data; % 'blue'
  color_cabezal = block.DialogPrm(10).Data; % 'blue'
  color_containers = block.DialogPrm(11).Data; % '#808080'(gris)  /'cyan'
  altura_viga_testera = block.DialogPrm(12).Data; % 5.000 m.
  % params.x0_mapa_containers = block.DialogPrm(13).Data; %  0.862 m.
  % params.y0_mapa_containers = block.DialogPrm(14).Data; % -9.650 m.

  fig = figure('Name', 'Simulacion', 'NumberTitle', 'off', 'Color', 'white', 'WindowState', 'maximized');
  setappdata(0, 'myFigHandle', fig);

  % Se obtiene la ruta del archivo SIMULACION.m, y se busca la imagen de fondo:
  rutaImagen = fullfile(fileparts(mfilename('fullpath')), 'Fondo_Simulacion.jpg');
  img = imread(rutaImagen);

  %% Escalado Ejes Imagen:
  img = fliplr(img); %Voltear la imagen horizontalmente, porque set(gca, 'YDir', 'normal') luego la da vuelta.
  img = rot90(img, 2); %Voltear la imagen verticalmente, porque set(gca, 'YDir', 'normal') luego la da vuelta.
  anchoImagen=1246; % píxeles
  alturaImagen=890; % píxeles
  TreintaMetrosX=355; % píxeles (medidos en imagen)
  CuarentaYCincoMetrosY=471.5; % píxeles (medidos en imagen)
  FactorEscalaX=30.000/TreintaMetrosX; % 84.50704225
  FactorEscalaY=45.000/CuarentaYCincoMetrosY; % 95.44008484
  NuevoRangoX = round(anchoImagen*FactorEscalaX); % round(105.2957746) = 105.296 [píxeles/1000 y m]
  NuevoRangoY = round(alturaImagen*FactorEscalaY); % round(84.9416755) = 84.942 [píxeles/1000 y m]
  OrigenX=47.1505;% [píxeles/1000 y m]
  OrigenY=12.0868;% [píxeles/1000 y m]

  imshow(img, 'InitialMagnification', 'fit', 'XData', [-OrigenX, NuevoRangoX-OrigenX], 'YData', [-OrigenY, NuevoRangoY-OrigenY]);
  axis on; % Mostrar los ejes
  % Invertir eje Y para que aumente hacia arriba, y cambiar porcentaje de Posición[%x,%y,%ancho,%alto] del cuadro de la figura:
  set(gca, 'Box', 'on', 'YDir', 'normal','Position', [0.07, 0.07, 0.9, 0.9]);
  % %% BUQUE:
  % color_buque = 'black';
  % patch('Faces', [1, 2, 3], 'Vertices', [params.x0_mapa_containers, params.y0_mapa_containers+H_c*7; ...
  %   params.x0_mapa_containers+W_c, params.y0_mapa_containers+H_c*7; params.x0_mapa_containers+W_c, ...
  %   -10.800], 'FaceColor', color_buque, 'EdgeColor', 'none'); % Triángulo izquierdo
  % rectangle('Position',[3.300, -10.800, W_c*17, 1.150],'FaceColor', color_buque,'LineWidth', 0.3); % Rectángulo fondo
  % patch('Faces', [1, 2, 3], 'Vertices', [params.x0_mapa_containers+W_c*18, -10.800; params.x0_mapa_containers+W_c*18, ...
  %       params.y0_mapa_containers+H_c*7; params.x0_mapa_containers+W_c*19, params.y0_mapa_containers+H_c*7], ...
  %       'FaceColor', color_buque, 'EdgeColor', 'none'); % Triángulo derecho
  % 
  % %% CONTAINERS DE LOS CAMIONES:
  % % Creamos los objetos de los containers de los camiones. 5 carriles.
  % PosY_Carriles=evalin('base', 'posYcc');   % Coordenadas de la esquina inferior izquierda de los
  % PosX_Carriles=evalin('base', 'posX_I_C'); % containers de los carriles.
  % mapa_containers_carriles=evalin('base', 'mapa_containers_carriles');
  % containers_camiones_handles = gobjects(1,5); % Objetos gráficos de los containers de camiones.
  % for i=1:5
  %     containers_camiones_handles(1,i) = rectangle('Position',[PosX_Carriles(i), PosY_Carriles, W_c, H_c],'FaceColor', ...
  %                                               color_containers,'LineWidth', 0.3, 'EdgeColor', 'white', 'LineWidth', 1);
  %     if mapa_containers_carriles(i) == 1
  %       containers_camiones_handles(1,i).Visible = 'on';
  %     else
  %       containers_camiones_handles(1,i).Visible = 'off';
  %     end
  % end
  % %% CONTAINERS DEL BARCO:
  % % Creamos los objetos de los containers del barco.
  % mapa_containers = evalin('base', 'mapa_containers'); % Tomamos el array que se encuentra en el Workspace.
  % mapa_containers = mapa_containers(end:-1:1,:); % Se realiza una permutación de las filas, para que coincidan ...
  %                                                % los índices de la matriz con la numeración de los containers.
  % [n_filas, n_columnas] = size(mapa_containers); % [14, 19]
  % containers_barco_handles = gobjects(n_columnas, n_filas); % Handles de los containers del barco.
  % 
  % for j=1:n_columnas
  %   for i=1:n_filas
  %       containers_barco_handles(j, i) = rectangle('Position',[params.x0_mapa_containers+W_c*(j-1), ...
  %                                                   params.y0_mapa_containers+H_c*(i-1), W_c, H_c], ...
  %                                                   'FaceColor', color_containers,'LineWidth', 0.3, ...
  %                                                   'EdgeColor', 'white', 'LineWidth', 1);
  %       if mapa_containers(i, j) == 1
  %         containers_barco_handles(j, i).Visible = 'on';
  %       else
  %         containers_barco_handles(j, i).Visible = 'off';
  %       end
  %   end
  % end
  % 
  % %% LÍNEAS TRANSICIÓN MODO MANUAL/AUTO:
  % lineas_horizontales_handles = gobjects(22); % Handles de las líneas de transición.
  % lineas_verticales_handles = gobjects(22);   %
  % 
  % % Línea Horizontal:
  % lineas_horizontales_handles(1) = ...
  %   line([-47.1928              params.x0_mapa_containers+W_c*(1/2-1)],  ...
  %        [ altura_viga_testera  altura_viga_testera]...
  %        ,'Color','#E99000','LineWidth',1);
  % 
  % max_cant_c_apilados = nnz(mapa_containers(:,1)==-1|mapa_containers(:,1)==1); %Comienza siendo siempre
  %                                                                     % la cant. de conts. de la 1er col.
  % % Línea Vertical:
  % lineas_verticales_handles(1) = ...
  %   line([params.x0_mapa_containers+W_c*(1/2-1)  params.x0_mapa_containers+W_c*(1/2-1)], ...
  %        [altura_viga_testera             params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)]...
  %        ,'Color','#E99000','LineWidth',1);
  % for i=1:n_columnas %[1, 19]
  %     % Líneas Horizontales:
  %     lineas_horizontales_handles(i+1) = ...
  %       line([params.x0_mapa_containers+W_c*(1/2+i-2)                  params.x0_mapa_containers+W_c*(1/2+i-1)], ...
  %            [params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)  params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)]...
  %            ,'Color','#E99000','LineWidth',1);
  %     aux_cant_c_apilados = max_cant_c_apilados;
  %     if(i+1 <= n_columnas)
  %       max_cant_c_apilados = max( nnz(mapa_containers(:,i)==-1|mapa_containers(:,i)==1), ...
  %                                  nnz(mapa_containers(:,i+1)==-1|mapa_containers(:,i+1)==1) );
  %       if( aux_cant_c_apilados > nnz(mapa_containers(:,i)==-1|mapa_containers(:,i)==1)  && ...
  %           aux_cant_c_apilados > nnz(mapa_containers(:,i+1)==-1|mapa_containers(:,i+1)==1) )
  %         if(i+2 <= n_columnas)
  %           if( aux_cant_c_apilados < nnz(mapa_containers(:,i+2)==-1|mapa_containers(:,i+2)==1) )
  %             max_cant_c_apilados = aux_cant_c_apilados;
  %           elseif( nnz(mapa_containers(:,i)==-1|mapa_containers(:,i)==1) < nnz(mapa_containers(:,i+2)==-1|mapa_containers(:,i+2)==1) )
  %             max_cant_c_apilados = nnz(mapa_containers(:,i+2)==-1|mapa_containers(:,i+2)==1);
  %           end
  %         end
  %       end
  %     end
  % 
  %     % Líneas Verticales:
  %     lineas_verticales_handles(i+1) = ...
  %       line([params.x0_mapa_containers+W_c*(1/2+i-1)                  params.x0_mapa_containers+W_c*(1/2+i-1)], ...
  %            [params.y0_mapa_containers+H_c*(aux_cant_c_apilados+1/2)  params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)]...
  %            ,'Color','#E99000','LineWidth',1);
  % end
  % 
  % % Línea Horizontal:
  % lineas_horizontales_handles(21) = ...
  %   line([params.x0_mapa_containers+W_c*(1/2+18)                   params.x0_mapa_containers+W_c*(1/2+19)], ...
  %        [params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)  params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)]...
  %        ,'Color','#E99000','LineWidth',1);
  % % Línea Vertical:
  % lineas_verticales_handles(21) = ...
  %   line([params.x0_mapa_containers+W_c*(1/2+19)                   params.x0_mapa_containers+W_c*(1/2+19)], ...
  %        [params.y0_mapa_containers+H_c*(max_cant_c_apilados+1/2)  -11.400]...
  %        ,'Color','#E99000','LineWidth',1);

  %% CARRO:
  carro_handle = rectangle('Position', [x_t, 45.000, ancho_carro, alto_carro], 'FaceColor', color_carro, ...
                          'EdgeColor', 'white', 'LineWidth', 1);
  %% CABLE:
  wirerope_handle = line([x_l x_t], [(y_l+alto_spreader+alto_headblock) 45000],'Color','black','LineWidth',1.5);
  %% CABEZAL:
  % Create hgtransform objects for rotation:
  hTransHeadblock = hgtransform('Parent', fig.CurrentAxes);
  hTransSpreader  = hgtransform('Parent', fig.CurrentAxes);
  hTransContainer = hgtransform('Parent', fig.CurrentAxes);

  headblock_handle = rectangle('Position',[(x_l-ancho_headblock/2),y_l+alto_spreader,ancho_headblock, ...
                              alto_headblock],'FaceColor', color_cabezal, 'EdgeColor', 'white', 'LineWidth', 0.01, ...
                              'Parent', hTransHeadblock);
  spreader_handle = rectangle('Position',[(x_l-ancho_spreader/2),y_l,ancho_spreader, ...
                              alto_spreader],'FaceColor', color_cabezal, 'EdgeColor', 'white', 'LineWidth', 0.01, ...
                              'Parent', hTransSpreader);
  container_handle = rectangle('Position',[x_l-W_c/2, y_l-H_c, W_c, H_c], ...
                              'FaceColor', color_containers,'LineWidth', 0.3, 'EdgeColor', 'white', 'LineWidth', 1, ...
                              'Parent', hTransContainer);
  container_handle.Visible = 'off'; % Empezamos siempre sin container.
 
  
  % %% TIEMPO:
  % text_clock_handle = text(60.000, 72.000, sprintf('Tiempo: %.1f', clock), ...
  %                          'FontSize', 10, 'Color', 'black', 'BackgroundColor', 'white');
  % %% TESTIGOS LUMINOSOS:
  % manualPos = [-46.4, -3.3, 22.5, 1.9]; %[x_inicial, y_inicial, ancho, alto]
  % manual_handle = rectangle('Position', manualPos,'FaceColor', 'green', 'EdgeColor', 'black', 'LineWidth', 1);
  % text_manual_handle = text( manualPos(1) + manualPos(3)/2, ...  % centrado en x
  %                            manualPos(2) + manualPos(4)/2, ...  % centrado en y
  %                            'Modo Manual', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
  %                            'FontSize',10, 'FontWeight','bold', 'Color','black' );
  % automaticoPos = [-23.9, -3.3, 22.5, 1.9];
  % automatico_handle = rectangle('Position', automaticoPos,'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'black', 'LineWidth', 1);
  % text_automatico_handle = text( automaticoPos(1) + automaticoPos(3)/2, ...
  %                                automaticoPos(2) + automaticoPos(4)/2, ...
  %                               'Modo Automático', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
  %                               'FontSize',10, 'FontWeight','bold', 'Color',[0.7 0.7 0.7] );
  % 
  % twistlocksAbiertosPos = [-46.4, -3.3-1.9*1, 22.5, 1.9]; %[x_inicial, y_inicial, ancho, alto]
  % twistlocks_abiertos_handle = rectangle('Position', twistlocksAbiertosPos,'FaceColor', 'green', 'EdgeColor', 'black', 'LineWidth', 1);
  % text_twistlocks_abiertos_handle = text( twistlocksAbiertosPos(1) + twistlocksAbiertosPos(3)/2, ...  % centrado en x
  %                                        twistlocksAbiertosPos(2) + twistlocksAbiertosPos(4)/2, ...  % centrado en y
  %                                        'Twistlocks Abiertos', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
  %                                        'FontSize',10, 'FontWeight','bold', 'Color','black' );
  % twistlocksCerradosPos = [-23.9, -3.3-1.9*1, 22.5, 1.9];
  % twistlocks_cerrados_handle = rectangle('Position', twistlocksCerradosPos,'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'black', 'LineWidth', 1);
  % text_twistlocks_cerrados_handle = text( twistlocksCerradosPos(1) + twistlocksCerradosPos(3)/2, ...
  %                                         twistlocksCerradosPos(2) + twistlocksCerradosPos(4)/2, ...
  %                                         'Twistlocks Cerrados', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
  %                                         'FontSize',10, 'FontWeight','bold', 'Color',[0.7 0.7 0.7] );  
  % 
  % balanceoActivadoPos = [-46.4, -3.3-1.9*2, 22.5, 1.9]; %[x_inicial, y_inicial, ancho, alto]
  % balanceo_activado_handle = rectangle('Position', balanceoActivadoPos,'FaceColor', 'green', 'EdgeColor', 'black', 'LineWidth', 1);
  % text_balanceo_activado_handle = text( balanceoActivadoPos(1) + balanceoActivadoPos(3)/2, ...  % centrado en x
  %                                       balanceoActivadoPos(2) + balanceoActivadoPos(4)/2, ...  % centrado en y
  %                                       'Control Balanceo Activado', 'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
  %                                       'FontSize',10, 'FontWeight','bold', 'Color','black' );
  % balanceoDesactivadoPos = [-23.9, -3.3-1.9*2, 22.5, 1.9];
  % balanceo_desactivado_handle = rectangle('Position', balanceoDesactivadoPos,'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'black', 'LineWidth', 1);
  % text_balanceo_desactivado_handle = text( balanceoDesactivadoPos(1) + balanceoDesactivadoPos(3)/2, ...
  %                                          balanceoDesactivadoPos(2) + balanceoDesactivadoPos(4)/2, ...
  %                                          'Control Balanceo Desactivado', 'HorizontalAlignment','center', 'VerticalAlignment', ...
  %                                          'middle', 'FontSize',10, 'FontWeight','bold', 'Color',[0.7 0.7 0.7] );
  % 
  % cargaNominalPos = [-46.4, -3.3-1.9*3, 22.5, 1.9]; %[x_inicial, y_inicial, ancho, alto]
  % carga_nominal_handle = rectangle('Position', cargaNominalPos,'FaceColor', 'green', 'EdgeColor', 'black', 'LineWidth', 1);
  % text_carga_nominal_handle = text( cargaNominalPos(1) + cargaNominalPos(3)/2, ...  % centrado en x
  %                                   cargaNominalPos(2) + cargaNominalPos(4)/2, ...  % centrado en y
  %                                   sprintf('Carga Nominal: %.0f kg', m_l_estimada), 'HorizontalAlignment','center', ... 
  %                                   'VerticalAlignment','middle', 'FontSize',10, 'FontWeight','bold', 'Color','black' );
  % SobrecargaPos = [-23.9, -3.3-1.9*3, 22.5, 1.9];
  % sobrecarga_handle = rectangle('Position', SobrecargaPos,'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'black', 'LineWidth', 1);
  % text_sobrecarga_handle = text( SobrecargaPos(1) + SobrecargaPos(3)/2, ...
  %                                SobrecargaPos(2) + SobrecargaPos(4)/2, ...
  %                                'Sobrecarga', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',10, ...
  %                                'FontWeight','bold', 'Color',[0.7 0.7 0.7] );
  % 
  % funcionamientoNormalPos = [-46.4, -3.3-1.9*4, 22.5, 1.9]; %[x_inicial, y_inicial, ancho, alto]
  % funcionamiento_normal_handle = rectangle('Position', funcionamientoNormalPos,'FaceColor', 'green', 'EdgeColor', 'black', 'LineWidth', 1);
  % text_funcionamiento_normal_handle = text( funcionamientoNormalPos(1) + funcionamientoNormalPos(3)/2, ...  % centrado en x
  %                                           funcionamientoNormalPos(2) + funcionamientoNormalPos(4)/2, ...  % centrado en y
  %                                           'Funcionamiento Normal', 'HorizontalAlignment','center', 'VerticalAlignment', ...
  %                                           'middle', 'FontSize',10, 'FontWeight','bold', 'Color','black' );
  % EmergenciaPos = [-23.9, -3.3-1.9*4, 22.5, 1.9];
  % emergencia_handle = rectangle('Position', EmergenciaPos,'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'black', 'LineWidth', 1);
  % text_emergencia_handle = text( EmergenciaPos(1) + EmergenciaPos(3)/2, ...
  %                                EmergenciaPos(2) + EmergenciaPos(4)/2, ...
  %                                'Emergencia', 'HorizontalAlignment','center', 'VerticalAlignment','middle', 'FontSize',10, ...
  %                                'FontWeight','bold', 'Color',[0.7 0.7 0.7] );


  % %% Pasaje de información entre funciones:
  fig.UserData = struct( ...
    'carro_handle', carro_handle, ...
    'wirerope_handle', wirerope_handle, ...
    'headblock_handle', headblock_handle, ...
    'spreader_handle', spreader_handle, ...
    'container_handle', container_handle, ...
    'hTransHeadblock', hTransHeadblock, ...
    'hTransSpreader',  hTransSpreader,  ...
    'hTransContainer', hTransContainer ...
);

  block.Dwork(1).Data = initialized; % Vector Dwork (Differential Work)
  % fig.UserData = struct('containers_barco_handles', containers_barco_handles, ...
  %                       'containers_camiones_handles', containers_camiones_handles, ...
  %                       'carro_handle', carro_handle, ...
  %                       'wirerope_handle', wirerope_handle, ...
  %                       'headblock_handle', headblock_handle, ...
  %                       'spreader_handle', spreader_handle, ...
  %                       'container_handle', container_handle, ...
  %                       'hTransHeadblock', hTransHeadblock, ...
  %                       'hTransSpreader', hTransSpreader, ...
  %                       'hTransContainer', hTransContainer, ...
  %                       'lineas_horizontales_handles', lineas_horizontales_handles, ...
  %                       'lineas_verticales_handles', lineas_verticales_handles, ...
  %                       'text_clock_handle', text_clock_handle, ...
  %                       'manual_handle', manual_handle, ...
  %                       'text_manual_handle', text_manual_handle, ...
  %                       'automatico_handle', automatico_handle, ...
  %                       'text_automatico_handle', text_automatico_handle, ...
  %                       'twistlocks_abiertos_handle', twistlocks_abiertos_handle,...
  %                       'text_twistlocks_abiertos_handle', text_twistlocks_abiertos_handle, ...
  %                       'twistlocks_cerrados_handle', twistlocks_cerrados_handle,...
  %                       'text_twistlocks_cerrados_handle', text_twistlocks_cerrados_handle, ...
  %                       'balanceo_activado_handle', balanceo_activado_handle,...
  %                       'text_balanceo_activado_handle', text_balanceo_activado_handle, ...
  %                       'balanceo_desactivado_handle', balanceo_desactivado_handle,...
  %                       'text_balanceo_desactivado_handle', text_balanceo_desactivado_handle, ...
  %                       'carga_nominal_handle', carga_nominal_handle,...
  %                       'text_carga_nominal_handle', text_carga_nominal_handle, ...
  %                       'sobrecarga_handle', sobrecarga_handle,...
  %                       'text_sobrecarga_handle', text_sobrecarga_handle, ...
  %                       'funcionamiento_normal_handle', funcionamiento_normal_handle,...
  %                       'text_funcionamiento_normal_handle', text_funcionamiento_normal_handle, ...
  %                       'emergencia_handle', emergencia_handle, ...
  %                       'text_emergencia_handle', text_emergencia_handle);
  % block.Dwork(1).Data = initialized; % Vector Dwork (Differential Work)
end %Start

% %%
% %% InitializeConditions:
% %%   Functionality    : Called at the start of simulation and if it is
% %%                      present in an enabled subsystem configured to reset
% %%                      states, it will be called when the enabled subsystem
% %%                      restarts execution to reset the states.
% %%   Required         : No
% %%   C MEX counterpart: mdlInitializeConditions
% function InitializeConditions(block)
%     %% Initialize Dwork:
%     % block.Dwork(1).Data = ;
% end %InitializeConditions

%%
%% Update:
%%   Functionality    : Called to update discrete states
%%                      during simulation step
%%   Required         : No
%%   C MEX counterpart: mdlUpdate
% function Update(block)
%
% end %Update

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C MEX counterpart: mdlOutputs
function Outputs(block)
  persistent fig params N_containers_anterior;

  x_t = block.InputPort(1).Data;
  x_l = block.InputPort(2).Data;
  y_l = block.InputPort(3).Data;
  theta_l = block.InputPort(4).Data; %[rad]
  % clock = block.InputPort(13).Data;
  % actuador_TWISLOCKS = block.InputPort(6).Data;
  % N_containers = block.InputPort(7).Data;
  % AlturasLineasTransicion = block.InputPort(8).Data;
  % MODO_MANUAL = block.InputPort(9).Data;
  % CONTROL_DE_BALANCEO = block.InputPort(10).Data;
  % CARGA_NOMINAL = block.InputPort(11).Data;
  % EMERGENCIA = block.InputPort(12).Data;
  % m_l_estimada = block.InputPort(5).Data;

  initialized = block.Dwork(1).Data;
  if (initialized == false) % Solo se entra la primera vez:
    %% Inicialización:
    initialized = true;
    block.Dwork(1).Data = initialized; % initialized = false.

    params.W_c = block.DialogPrm(1).Data; % Ancho de conteiner (estándar): 2.438 m.
    params.H_c = block.DialogPrm(2).Data; % Alto de conteiner (estándar): 2.591 m.
    params.alto_carro = block.DialogPrm(3).Data; % 1.500 m.
    params.alto_headblock = block.DialogPrm(4).Data; % 0.850 m.
    params.alto_spreader = block.DialogPrm(5).Data; % 0.750 m.
    params.ancho_carro = block.DialogPrm(6).Data; % 7.500 m.
    params.ancho_headblock = block.DialogPrm(7).Data; % 2.100 m.
    params.ancho_spreader = block.DialogPrm(8).Data; % 25.00 m.
    % params.x0_mapa_containers = block.DialogPrm(13).Data; %  0.862 m.
    % params.y0_mapa_containers = block.DialogPrm(14).Data; % -9.650 m.

    fig = getappdata(0, 'myFigHandle'); % Retrieve the figure handle

    %N_containers_anterior = N_containers; %zeros(1,24);
  end

  %% ACTUALIZACIÓN DE LOS DIBUJOS MÓVILES:
  if ~isempty(fig) && ishandle(fig) %Si se ha creado al figura y no se la ha cerrado:
    ud = fig.UserData;  % Cachea UserData
    %% CARRO:
    set(ud.carro_handle, 'Position', [x_t-(3/4)*params.ancho_carro, 45.000, params.ancho_carro, params.alto_carro]);

    %% CABEZAL:
    % Update rectangle positions (nonrotated positions):
    headblock_pos = [(x_l-params.ancho_headblock/2), y_l+params.alto_spreader, params.ancho_headblock, params.alto_headblock];
    spreader_pos = [(x_l-params.ancho_spreader/2), y_l, params.ancho_spreader, params.alto_spreader];
    set(ud.headblock_handle, 'Position', headblock_pos);
    set(ud.spreader_handle, 'Position', spreader_pos);

    % Compute the center of rotation [Xcentro,Ycentro,Zcentro], for each rectangle:
    center_headblock = [x_l, spreader_pos(2), 0];
    center_spreader  = [x_l, spreader_pos(2), 0];

    % Create transformation matrices to rotate about the center by theta_l (in radians):
    T_head = makehgtform('translate', center_headblock) * makehgtform('zrotate', theta_l) * ...
             makehgtform('translate', -center_headblock);
    T_spreader = makehgtform('translate', center_spreader) * makehgtform('zrotate', theta_l) * ...
                 makehgtform('translate', -center_spreader);

    % Apply the transformation matrices to the corresponding hgtransform objects:
    set(ud.hTransHeadblock, 'Matrix', T_head);
    set(ud.hTransSpreader, 'Matrix', T_spreader);

    % %% CONTAINER AGARRADO:
    % if actuador_TWISLOCKS
    %     container_pos = [(x_l-params.W_c/2), y_l-params.H_c, params.W_c, params.H_c];
    %     set(ud.container_handle, 'Position', container_pos);
    %     center_container = [x_l, spreader_pos(2),0];
    %     T_container = makehgtform('translate', center_container) * makehgtform('zrotate', theta_l) * ...
    %                   makehgtform('translate', -center_container);
    %     set(ud.hTransContainer, 'Matrix', T_container);
    %     ud.container_handle.Visible = 'on';
    % else
    %    ud.container_handle.Visible = 'off';
    % end
    %% CABLE:
    set(ud.wirerope_handle, 'XData', [x_l x_t], 'YData', [y_l 45.000]);

    % %% TIEMPO:
    % set(ud.text_clock_handle, 'String', sprintf('Tiempo: %.1f s', clock));
    % %% TESTIGOS LUMINOSOS:
    % if MODO_MANUAL
    %     set(ud.manual_handle, 'FaceColor','green');
    %     set(ud.text_manual_handle, 'Color','black');
    %     set(ud.automatico_handle,'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_automatico_handle,'Color',[0.7 0.7 0.7])
    % else
    %     set(ud.manual_handle, 'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_manual_handle, 'Color',[0.7 0.7 0.7]);
    %     set(ud.automatico_handle, 'FaceColor','green');
    %     set(ud.text_automatico_handle, 'Color','black');
    % end
    % 
    % if ~actuador_TWISLOCKS
    %     set(ud.twistlocks_abiertos_handle, 'FaceColor','green');
    %     set(ud.text_twistlocks_abiertos_handle, 'Color','black');
    %     set(ud.twistlocks_cerrados_handle,'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_twistlocks_cerrados_handle,'Color',[0.7 0.7 0.7]);
    % else
    %     set(ud.twistlocks_abiertos_handle, 'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_twistlocks_abiertos_handle, 'Color',[0.7 0.7 0.7]);
    %     set(ud.twistlocks_cerrados_handle, 'FaceColor','green');
    %     set(ud.text_twistlocks_cerrados_handle, 'Color','black');
    % end
    % 
    % if CONTROL_DE_BALANCEO
    %     set(ud.balanceo_activado_handle, 'FaceColor','green');
    %     set(ud.text_balanceo_activado_handle, 'Color','black');
    %     set(ud.balanceo_desactivado_handle,'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_balanceo_desactivado_handle,'Color',[0.7 0.7 0.7]);
    % else
    %     set(ud.balanceo_activado_handle, 'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_balanceo_activado_handle, 'Color',[0.7 0.7 0.7]);
    %     set(ud.balanceo_desactivado_handle, 'FaceColor','green');
    %     set(ud.text_balanceo_desactivado_handle, 'Color','black');
    % end
    % 
    % if CARGA_NOMINAL
    %     set(ud.carga_nominal_handle, 'FaceColor','green');
    %     set(ud.text_carga_nominal_handle, 'Color','black', 'String',sprintf('Carga Nominal: %.0f kg', m_l_estimada));
    %     set(ud.sobrecarga_handle,'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_sobrecarga_handle,'Color',[0.7 0.7 0.7], 'String','Sobrecarga');
    % else
    %     set(ud.carga_nominal_handle, 'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_carga_nominal_handle, 'Color',[0.7 0.7 0.7], 'String','Carga Nominal');
    %     set(ud.sobrecarga_handle, 'FaceColor','#E99000');
    %     set(ud.text_sobrecarga_handle, 'Color','black', 'String',sprintf('Sobrecarga: %.0f kg', m_l_estimada));
    % end
    % 
    % if ~EMERGENCIA
    %     set(ud.funcionamiento_normal_handle, 'FaceColor','green');
    %     set(ud.text_funcionamiento_normal_handle, 'Color','black');
    %     set(ud.emergencia_handle,'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_emergencia_handle,'Color',[0.7 0.7 0.7]);
    % else
    %     set(ud.funcionamiento_normal_handle, 'FaceColor',[0.9 0.9 0.9]);
    %     set(ud.text_funcionamiento_normal_handle, 'Color',[0.7 0.7 0.7]);
    %     set(ud.emergencia_handle, 'FaceColor','red');
    %     set(ud.text_emergencia_handle, 'Color','black');
    % end

    %% 
    % drawnow limitrate; % drawnow
    % 
    % if ~isequal(N_containers, N_containers_anterior) % Solo si se modifica algún container:
    %   % changed_index = find(N_containers ~= N_containers_anterior);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    %   %% CONTAINERS EN REPOSO:
    %   % En los Carriles de los Camiones:
    %   for i=1:5
    %     if(N_containers(i) > 0)
    %       ud.containers_camiones_handles(1,i).Visible = 'on';
    %     else
    %       ud.containers_camiones_handles(1,i).Visible = 'off';
    %     end
    %   end
    %   % En los Bordes (Paredes) del Barco:
    %   for i=1:7
    %       % Lado izquierdo.
    %       if (i <= N_containers(6))
    %         ud.containers_barco_handles(1, i+7).Visible = 'on';
    %       else
    %         ud.containers_barco_handles(1, i+7).Visible = 'off';
    %       end
    %       % Lado derecho
    %       if (i <= N_containers(24))
    %         ud.containers_barco_handles(19, i+7).Visible = 'on';
    %       else
    %         ud.containers_barco_handles(19, i+7).Visible = 'off';
    %       end
    %   end
    %   % En la Parte Central del Barco:
    %   for j=2:18
    %     for i=1:14
    %       if (i <= N_containers(5+j))
    %         ud.containers_barco_handles(j, i).Visible = 'on';
    %       else
    %         ud.containers_barco_handles(j, i).Visible = 'off';
    %       end
    %     end
    %   end
    % 
    %   %% LÍNEAS TRANSICIÓN MODO MANUAL/AUTO:
    %   for i=1:(length(AlturasLineasTransicion)-1)
    %       % Líneas Horizontales:
    %       ud.lineas_horizontales_handles(i).YData = [AlturasLineasTransicion(i)
    %                                                  AlturasLineasTransicion(i)];
    %       % Líneas Verticales:
    %       ud.lineas_verticales_handles(i).YData = [AlturasLineasTransicion(i)
    %                                                  AlturasLineasTransicion(i+1)];
    %   end
    % end
    % N_containers_anterior = N_containers;
  else
    return;
  end

end %Outputs

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C MEX counterpart: mdlDerivatives
% function Derivatives(block)
%
% end %Derivatives

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C MEX counterpart: mdlTerminate
function Terminate(~)
  % block.Dwork(1).Data = false; % initialized = false.
end %Terminate






