clearvars;
clc;
%% ---------------------------------------------------- CREATE 'vem' STRUCT
thickness  = 1;
numElement = 5000;
[Node,Element,Supp,Load] = PolyMesher(@HookDomain,numElement,30);

Ae = zeros(numElement,1);
for iel = 1:numElement
    enodes  =   Element{iel};
    verts   =   Node(enodes,:);
    parea   =   polyarea(verts(:,1),verts(:,2));
    Ae(iel) =   parea;
end

vem = struct(...
  'NNode',size(Node,1),...     % Number of nodes
  'NElem',size(Element,1),...  % Number of elements
  'Node',Node,...              % [NNode x 2] array of nodes
  'Element',{Element},...      % [NElement x Var] cell array of elements
  'AreaElem',Ae,...            % Array of element area 
  'Thickness',thickness,...    % Thickness of the domain
  'Supp',Supp,...              % Array of supports
  'Load',Load,...              % Array of loads
  'Nu0',0.3,...                % Poisson's ratio of solid material
  'E0',1.0...                  % Young's modulus of solid material
   );
%% ---------------------------------------------------- CREATE 'opt' STRUCT
R = 2; P = Filter(Element, Node, R);
VolFrac = 0.4;
xIni = VolFrac*ones(size(P,2),1);
opt = struct(...               
  'xMin',0.0,...               % Lower bound for design variables
  'xMax',1.0,...               % Upper bound for design variables
  'xIni',xIni,...              % Initial design variables
  'P',P,...                    % Matrix that uses as filter.
  'VolFrac',VolFrac,...        % Specified volume fraction cosntraint
  'VolElem', Ae*thickness,...  % Array of element volume 
  'Tol',0.01,...               % Convergence tolerance on design vars.
  'MaxIter',100,...            % Max. number of optimization iterations
  'Alpha',0.6...               % History parameter
   );              
%% ---------------------------------------------------------- RUN 'VTop'
figure;
for penal = 1:0.5:4            % Continuation on the penalty parameter
   disp(['current p: ', num2str(penal)]);
   opt.MatIntFnc = @(y)MatIntFnc(y,'SIMP',penal); % Material interpolation.
   [opt.xIni,vem] = VTop(vem,opt);
end
%% ------------------------------------------------------------------------