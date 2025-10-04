function Simulacion(block)
    setup(block);
end

%% ========================== setup =================================
function setup(block)
    block.NumInputPorts  = 5;
    block.NumOutputPorts = 0;

    block.SetPreCompInpPortInfoToDynamic;

    % x_t
    block.InputPort(1).Dimensions        = 1;
    block.InputPort(1).DatatypeID        = 0;
    block.InputPort(1).Complexity        = 'Real';
    block.InputPort(1).SamplingMode      = 'Sample';
    block.InputPort(1).DirectFeedthrough = true;

    % x_l
    block.InputPort(2).Dimensions        = 1;
    block.InputPort(2).DatatypeID        = 0;
    block.InputPort(2).Complexity        = 'Real';
    block.InputPort(2).SamplingMode      = 'Sample';
    block.InputPort(2).DirectFeedthrough = true;

    % y_l
    block.InputPort(3).Dimensions        = 1;
    block.InputPort(3).DatatypeID        = 0;
    block.InputPort(3).Complexity        = 'Real';
    block.InputPort(3).SamplingMode      = 'Sample';
    block.InputPort(3).DirectFeedthrough = true;

    % theta_l
    block.InputPort(4).Dimensions        = 1;
    block.InputPort(4).DatatypeID        = 0;
    block.InputPort(4).Complexity        = 'Real';
    block.InputPort(4).SamplingMode      = 'Sample';
    block.InputPort(4).DirectFeedthrough = true;

    % Perfil de obstáculos (vector 1-D de 33)
    block.InputPort(5).Dimensions        = 33;
    block.InputPort(5).DatatypeID        = 0;
    block.InputPort(5).Complexity        = 'Real';
    block.InputPort(5).SamplingMode      = 'Sample';
    block.InputPort(5).DirectFeedthrough = true;

    block.NumDialogPrms     = 12;
    block.DialogPrmsTunable = repmat({'Nontunable'},1,block.NumDialogPrms);

    block.SampleTimes = [0.1 0];
    block.SimStateCompliance = 'DefaultSimState';

    block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
    block.RegBlockMethod('Start',                @Start);
    block.RegBlockMethod('Outputs',              @Outputs);
    block.RegBlockMethod('Terminate',            @Terminate);
end

%% ==================== DoPostPropSetup =============================
function DoPostPropSetup(block)
    block.NumDworks                 = 1;
    block.Dwork(1).Name            = 'initialized';
    block.Dwork(1).Dimensions      = 1;
    block.Dwork(1).DatatypeID      = 8;    % boolean
    block.Dwork(1).Complexity      = 'Real';
    block.Dwork(1).UsedAsDiscState = true;
end

%% ============================ Start ================================
function Start(block)
  close all;
  persistent fig img;

  initialized = false;

  x_t = block.InputPort(1).Data;
  x_l = block.InputPort(2).Data;
  y_l = block.InputPort(3).Data;

  W_c                 = block.DialogPrm(1).Data;
  H_c                 = block.DialogPrm(2).Data;
  alto_carro          = block.DialogPrm(3).Data;
  alto_headblock      = block.DialogPrm(4).Data;
  alto_spreader       = block.DialogPrm(5).Data;
  ancho_carro         = block.DialogPrm(6).Data;
  ancho_headblock     = block.DialogPrm(7).Data;
  ancho_spreader      = block.DialogPrm(8).Data;
  color_carro         = block.DialogPrm(9).Data;
  color_cabezal       = block.DialogPrm(10).Data;
  color_containers    = block.DialogPrm(11).Data;
  altura_viga_testera = block.DialogPrm(12).Data; %#ok<NASGU>

  fig = figure('Name','Simulacion','NumberTitle','off','Color','white','WindowState','maximized');
  setappdata(0,'myFigHandle',fig);

  rutaImagen = fullfile(fileparts(mfilename('fullpath')), 'Fondo_Simulacion.jpg');
  img = imread(rutaImagen);

  img = fliplr(img);
  img = rot90(img,2);
  anchoImagen           = 1246;
  alturaImagen          = 890;
  TreintaMetrosX        = 355;
  CuarentaYCincoMetrosY = 471.5;
  FactorEscalaX         = 30.000/TreintaMetrosX;
  FactorEscalaY         = 45.000/CuarentaYCincoMetrosY;
  NuevoRangoX           = round(anchoImagen*FactorEscalaX);
  NuevoRangoY           = round(alturaImagen*FactorEscalaY);
  OrigenX               = 47.1505;
  OrigenY               = 12.0868;

  imshow(img,'InitialMagnification','fit','XData',[-OrigenX, NuevoRangoX-OrigenX],'YData',[-OrigenY, NuevoRangoY-OrigenY]);
  axis on;
  set(gca,'Box','on','YDir','normal','Position',[0.07, 0.07, 0.9, 0.9]);

  % Carro
  carro_handle = rectangle('Position',[x_t, 45.000, ancho_carro, alto_carro], ...
                           'FaceColor',color_carro,'EdgeColor','white','LineWidth',1);
  % Cable
  wirerope_handle = line([x_l x_t], [(y_l+alto_spreader+alto_headblock) 45000], 'Color','black','LineWidth',1.5);

  % Cabezal con transform
  hTransHeadblock = hgtransform('Parent', fig.CurrentAxes);
  hTransSpreader  = hgtransform('Parent', fig.CurrentAxes);
  hTransContainer = hgtransform('Parent', fig.CurrentAxes);

  headblock_handle = rectangle('Position',[(x_l-ancho_headblock/2), y_l+alto_spreader, ancho_headblock, alto_headblock], ...
                               'FaceColor',color_cabezal,'EdgeColor','white','LineWidth',0.01,'Parent',hTransHeadblock);
  spreader_handle  = rectangle('Position',[(x_l-ancho_spreader/2),  y_l,                 ancho_spreader, alto_spreader], ...
                               'FaceColor',color_cabezal,'EdgeColor','white','LineWidth',0.01,'Parent',hTransSpreader);
  container_handle = rectangle('Position',[x_l-W_c/2, y_l-H_c, W_c, H_c], ...
                               'FaceColor',color_containers,'EdgeColor','white','LineWidth',1,'Parent',hTransContainer);
  container_handle.Visible = 'off';

  % Perfil de obstáculos (señal cuadrada con stairs)
  ax = fig.CurrentAxes;
  xlims = get(ax,'XLim');

  y0 = block.InputPort(5).Data;             % puede venir vacío al inicio
  if isempty(y0)
      NOBS = 33;
      y0 = zeros(1,NOBS);
  else
      y0  = double(y0(:))';
      NOBS = numel(y0);
  end
  x_centros = linspace(xlims(1)+2, xlims(2)-2, NOBS);
  [Xs, Ys]  = stairs(x_centros, y0);
  hold(ax,'on');
  obstacles_handle = plot(ax, Xs, Ys, 'LineWidth',2,'Color','#E99000');

  fig.UserData = struct( ...
      'carro_handle',     carro_handle, ...
      'wirerope_handle',  wirerope_handle, ...
      'headblock_handle', headblock_handle, ...
      'spreader_handle',  spreader_handle, ...
      'container_handle', container_handle, ...
      'hTransHeadblock',  hTransHeadblock, ...
      'hTransSpreader',   hTransSpreader, ...
      'hTransContainer',  hTransContainer, ...
      'obstacles_handle', obstacles_handle, ...
      'x_centros',        x_centros ...
  );

  block.Dwork(1).Data = initialized;
end

%% =========================== Outputs ===============================
function Outputs(block)
  persistent fig params

  x_t          = block.InputPort(1).Data;
  x_l          = block.InputPort(2).Data;
  y_l          = block.InputPort(3).Data;
  theta_l      = block.InputPort(4).Data;
  alturas_obst = double(block.InputPort(5).Data(:))';   % fila

  initialized = block.Dwork(1).Data;
  if (initialized == false)
      initialized = true;
      block.Dwork(1).Data = initialized;

      params.W_c             = block.DialogPrm(1).Data;
      params.H_c             = block.DialogPrm(2).Data;
      params.alto_carro      = block.DialogPrm(3).Data;
      params.alto_headblock  = block.DialogPrm(4).Data;
      params.alto_spreader   = block.DialogPrm(5).Data;
      params.ancho_carro     = block.DialogPrm(6).Data;
      params.ancho_headblock = block.DialogPrm(7).Data;
      params.ancho_spreader  = block.DialogPrm(8).Data;

      fig = getappdata(0,'myFigHandle');
  end

  if isempty(fig) || ~ishandle(fig), return; end
  ud = fig.UserData;
  if ~isstruct(ud), return; end

  % Carro
  if isfield(ud,'carro_handle')
      set(ud.carro_handle,'Position',[x_t-(3/4)*params.ancho_carro, 45.000, params.ancho_carro, params.alto_carro]);
  end

  % Cabezal (posiciones sin rotación)
  headblock_pos = [(x_l-params.ancho_headblock/2), y_l+params.alto_spreader, params.ancho_headblock, params.alto_headblock];
  spreader_pos  = [(x_l-params.ancho_spreader/2),  y_l,                      params.ancho_spreader, params.alto_spreader];
  if isfield(ud,'headblock_handle'), set(ud.headblock_handle,'Position',headblock_pos); end
  if isfield(ud,'spreader_handle'),  set(ud.spreader_handle, 'Position',spreader_pos);  end

  % Rotación
  center = [x_l, spreader_pos(2), 0];
  T = makehgtform('translate',center) * makehgtform('zrotate',theta_l) * makehgtform('translate',-center);
  if isfield(ud,'hTransHeadblock'), set(ud.hTransHeadblock,'Matrix',T); end
  if isfield(ud,'hTransSpreader'),  set(ud.hTransSpreader, 'Matrix',T); end

  % Cable
  if isfield(ud,'wirerope_handle')
      set(ud.wirerope_handle,'XData',[x_l x_t],'YData',[y_l 45.000]);
  end

  % Perfil de obstáculos (stairs)
  if isfield(ud,'obstacles_handle')
      xc = ud.x_centros;
      if numel(alturas_obst) ~= numel(xc)
          ax    = fig.CurrentAxes;
          xlims = get(ax,'XLim');
          xc    = linspace(xlims(1)+2, xlims(2)-2, numel(alturas_obst));
          ud.x_centros = xc;
          set(fig,'UserData',ud);
      end
      [Xs, Ys] = stairs(xc, alturas_obst);
      set(ud.obstacles_handle, 'XData', Xs, 'YData', Ys);
  end

  drawnow limitrate;
end

%% ========================== Terminate =============================
function Terminate(~)
end
