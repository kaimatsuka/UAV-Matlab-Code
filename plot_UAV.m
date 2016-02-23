function [h] = plot_UAV(wing,htail,vtail,fuse,prop) 
%C130 draws a simple 3D airplane loosely modelled after a Lockheed C-130.
% 
%% c130 Syntax 
% 
%  c130 
%  c130(x,y,z)
%  c130(...,'roll',RollDegrees)
%  c130(...,'pitch',PitchDegrees)
%  c130(...,'yaw',YawDegrees)
%  c130(...,'color',AirplaneColor)
%  c130(...,'fuselage',FuseLageColor)
%  c130(...,'wing',WingColor)
%  c130(...,'tailwing',TailwingColor)
%  c130(...,'fin',FinColor)
%  c130(...,'prop',PropellerColor)
%  c130(...,'scale',SizeScaleFactor)
%  c130(...,'z',ZScaleFactor)
%  c130(...,'linestyle','LineStyle')
%  c130(...,'linecolor',LineColor)
%  h = c130(...)
% 
%% c130 Description
%
% c130 draws a 3D airplane.
% 
% c130(x,y,z) draws an airplane centered approximately at the location
% given by x,y,z, where x, y, and z must be scalar values. 
% 
% c130(...,'roll',RollDegrees) specifies roll of the aircraft in degrees
% relative to the approximate center of gravity of the aircraft. 
%
% c130(...,'pitch',PitchDegrees) specifies pitch of the aircraft in
% degrees.
%
% c130(...,'yaw',YawDegrees) specifies yaw of the aircraft in degrees in
% the global coordinate system. The xyz2rpy function may help with
% determining appropriate yaw values. 
%
% c130(...,'color',AirplaneColor) specifies a color of all surfaces of
% the plane.  Color may be given by Matlab color name (e.g., 'red'),
% Matlab color abbreviation (e.g., 'r'), or RGB value (e.g., [1 0 0]).
% The 'color' option may be paired with options for specifying colors of
% specific parts of the plane.  For example, c130('color','b','wing',y)
% creates a blue airplane with yellow wings. Default color is gray. 
%
% c130(...,'fuselage',FuselageColor) specifies fuselage color. 
%
% c130(...,'wing',WingColor) specifies wing color. 
%
% c130(...,'tailwing',TailwingColor) specifies color of horizontal
% stabilizer wings at the tail of the plane. 
%
% c130(...,'fin',FinColor) specifies color of the vertical stabilizer fin
% at the tail of the plane. 
%
% c130(...,'prop',PropellerColor) specifies propeller color. 
%
% c130(...,'scale',SizeScaleFactor) scales dimensions of the plane by
% some scalar factor. By default, c130 draws an airplane with dimensions
% which approximately match the dimensions C-130 airplane in meters.  If you'd 
% like to draw a C-130 in dimensions of feet, try c130('scale',3.281). 
%
% c130(...,'z',ZScaleFactor) scales only vertical dimensions by some
% scalar value.  This may be useful if you're animating an airplane in a
% coordinate system where the vertical dimension is stretched to show
% relief.  If your axes are not equal, consider a ZScaleFactor given by
% the ratio of (xlim(2)-xlim(1))/(zlim(2)-zlim(1)). 
%
% c130(...,'linestyle','LineStyle') specifies style of the lines
% connecting surface vertices. Default is '-', but can be set to any
% valid linestyle, including 'none' to show a smooth surface. 
%
% c130(...,'linecolor',LineColor) specifies edge color of lines
% connecting surface verticies. Default line color is black. 
% 
% h = c130(...) returns handles of the 13 surface objects created by
% c130. Handles correspond to the following: 
% 
%   * h(1): main fuselage. 
%   * h(2): nose. 
%   * h(3): tail section of fuselage. 
%   * h(4): top side of main wing. 
%   * h(5): under side of main wing. 
%   * h(6): top side of tail wing. 
%   * h(7): under side of tail wing. 
%   * h(8): right side of fin. 
%   * h(9): left side of fin. 
%   * h(10:13): propellers.
%
%% Author Info
% Chad A. Greene of the University of Texas Institute for Geophysics (UTIG)
% wrote this on Saturday, September 27, 2014. Feel free to visit Chad on
% over at http://www.chadagreene.com. 
% 
% See also xyz2rpy

%% Set Defaults: 

wingColor     = .5*[1 1 1];
fusColor      = .5*[1 1 1];
tailWingColor = .5*[1 1 1];
finColor      = .5*[1 1 1];
propColor     = .5*[1 1 1];
edgeColor     = 'k';
linestyle     = '-'; 
scale         = 1; 
zscale        = 1; 
setroll       = false; 
setpitch      = false; 
setyaw        = false; 
x             = 0;
y             = 0;
z             = 0;

%% Parse Inputs: 
%{
if nargin>2 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
    x = varargin{1}; 
    y = varargin{2}; 
    z = varargin{3}; 
    assert(isscalar(x)==1,'c130 input x must be a scalar.')
    assert(isscalar(y)==1,'c130 input y must be a scalar.')
    assert(isscalar(z)==1,'c130 input z must be a scalar.')
end

tmp = strncmpi(varargin,'col',3); 
if any(tmp)
    wingColor     = varargin{find(tmp)+1}; 
    fusColor      = varargin{find(tmp)+1}; 
    tailWingColor = varargin{find(tmp)+1}; 
    finColor      = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'fus',3)|strcmpi(varargin,'body'); 
if any(tmp)
    fusColor = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'wing',4); 
if any(tmp)
    wingColor = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'tail',4); 
if any(tmp)
    tailWingColor = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'fin',3); 
if any(tmp)
    finColor = varargin{find(tmp)+1}; 
end

tmp = strncmpi(varargin,'prop',4); 
if any(tmp)
    propColor = varargin{find(tmp)+1}; 
end

tmp = strcmpi(varargin,'scale');  
if any(tmp) 
    scale = varargin{find(tmp)+1}; 
    assert(isscalar(scale)==1,'It may seem redundant, but the scale must be a scalar.')
end

tmp = strncmpi(varargin,'z',1); 
if any(tmp)
    zscale = varargin{find(tmp)+1}; 
    assert(isscalar(zscale)==1,'It may seem redundant, but the z scale must be a scalar.')
end

tmp = strncmpi(varargin,'lines',5); 
if any(tmp)
    linestyle = varargin{find(tmp)+1}; 
    assert(isnumeric(linestyle)==0,'LineStyle cannot be numeric.')
end

tmp = strncmpi(varargin,'linecol',7)|strncmpi(varargin,'edgecol',3); 
if any(tmp)
    edgeColor = varargin{find(tmp)+1}; 
end

tmp = strcmpi(varargin,'roll')|strcmpi(varargin,'biscuits'); 
if any(tmp)
    setroll = true; 
    roll = varargin{find(tmp)+1}; 
    assert(isscalar(roll)==1,'Roll must be a scalar value in degrees.')
end

tmp = strcmpi(varargin,'pitch'); 
if any(tmp)
    setpitch = true; 
    pitch = varargin{find(tmp)+1}; 
    assert(isscalar(pitch)==1,'Pitch must be a scalar value in degrees.')
end

tmp = strcmpi(varargin,'yaw')|strcmpi(varargin,'yall'); % to accommodate our southern friends 
if any(tmp)
    setyaw = true; 
    yaw = varargin{find(tmp)+1}; 
    assert(isscalar(yaw)==1,'Yaw must be a scalar value in degrees.')
end
%}

%% Dimensions: 

% TODO LIST:
%   adjust height of wing, tails,
%   addjust max thickness line of wing, tails(optional)
%   propeller hub
%   work on vtail shape
%   work on nose shape
%   work on htail shape
%   create cone part of propeller

% Input List
%{
% fuse.L  = 4.2857; % baseUAV
fuse.L  = 4.6315;
fuse.D  = 1.0833; % fuselage depth (also diameter)

wing.b  = 10;         % wing span
wing.c_r = 1.0256;    % wing root chord
wing.c_t = 0.5128;    % wing tip chord
wing.h_l_r = 1.0296;  % wing root leading-edge location from nose
wing.h_l_t = 1.1578;  % wing tail leading-edge location from nose
wing.t_r = 0.1231;    % wing thickness at root
wing.t_t = 0.0615;    % wing thickness at tip

wing.c  = 0.7692; % average chord length
% ^^ get rid of this

htail.b     = 3;       % horizontal tail span
htail.c_r   = 1.0325;  % horizontal tail root chord
htail.c_t   = 0.5059;  % horizontal tail tip chord
htail.h_l_r = 3.5990;  % horizontal tail root leading-edge location from nose
htail.h_l_t = 3.7306;  % horizontal tail tail leading-edge location from nose
htail.t_r   = 0.1239;  % horizontal tail thickness at root
htail.t_t   = 0.0607;  % horizontal tail thickness at tip

vtail.b     = 1.5;    % vertical tail span
vtail.c_r   = 0.9616; % vertical tail root chord
vtail.c_t   = 0.5769; % vertical tail tip chord
vtail.h_l_r = 3.6167;  % horizontal tail root leading-edge location from nose
vtail.h_l_t = 3.7129;  % horizontal tail tail leading-edge location from nose
vtail.t_r   = 0.1154;  % horizontal tail thickness at root
vtail.t_t   = 0.0692;  % horizontal tail thickness at tip

prop.h = 4.6315; % propeller location
prop.D = 1.62; % propeller diameter (ft)

%}


fusLength    = (fuse.L-fuse.D/2)/2;  % only middle part of fuselage
tailLength   = fusLength; 
fusRadius    = fuse.D/2;  
chordlength  = wing.c; % changing chordlength will affect several dimensions


%{

scale      = 0.25; % size the aircraft to be approximately the same as Disaster Relief Drone

wingspan     = 40*scale; 
tailwingspan = 16*scale; 
tailLength   = 12.5*scale; 
fusLength    = 13*scale;  
fusRadius    = 2.2*scale;  
wingWidth    = 5*scale; % changing wingWidth will affect several dimensions
%}

% Center the dimensions: 
% y = y+9*wingWidth/5; 
% z = z+fusRadius; 

%% Determine if a figure is already open: 

initialHoldState = 0; 
SetView = isempty(get(gcf,'CurrentAxes'));
    
if ~SetView
    if ishold
        initialHoldState = 1; % documents initial hold state to reset it later 
    else 
        cla % if hold is on initially, clear the axes and make a new plot
        SetView = true;
    end      
end

%% Draw surfaces: 

% Fuselage: 
[xcf,zcf,ycf] = cylinder(fusRadius,40); 
h(1) = surface(xcf+x,y+ycf*fusLength+fuse.D/2,z+zcf*zscale,...
    'facecolor',fusColor,'linestyle',linestyle,...
    'edgecolor',edgeColor);
if ~initialHoldState
    hold on
end

% Nose: 
[xcn,zcn,ycn] = cylinder(fusRadius.*([1 .95 .9 .8 .7 .5 .5 .5 .5]).*(cos(linspace(0,pi/2,9)).^.2)); 
zcn(6:end,:) = zcn(6:end,:)-fusRadius/5; % lower the tip of nose down
% ycn = -ycn.*.7*chordlength; 
ycn = -ycn*fusRadius+fuse.D/2;
h(2) = surface(x+xcn,y+ycn,z+zcn*zscale,...
    'facecolor',fusColor,'linestyle',linestyle,...
    'edgecolor',edgeColor);

% Tail (aft half of fuselage):
x1 = xcf(1,:); 
x2 = .8*x1;% zeros(size(x1)); 
y1 = fusLength*ones(size(x1))+fuse.D/2; 
y2 = y1+tailLength; 
z1 = zcf(1,:); 
z2 = fusRadius*ones(size(z1)); 
h(3) = surface(x+[x1;x2],y+[y1;y2],z+[z1;z2]*zscale,...
    'facecolor',fusColor,'linestyle',linestyle,...
    'edgecolor',edgeColor);

% Wings: 
xw1 = -linspace(-wing.b/2,wing.b/2,20); % wing width
% yw1 = 1.4*chordlength + abs(xw1)/100;       % 
yw1 = wing.h_l_r+(wing.h_l_t-wing.h_l_r)*abs(linspace(-1,1,length(xw1))); % leading edge
% yw2 = zeros(size(xw1))+1.4*chordlength + chordlength/3; 
% yw3 = yw1 + chordlength-abs(xw1)/20;        % aft line
yw2 = wing.h_l_r+wing.c_r*0.3+zeros(size(xw1));
yw3 = yw1+wing.c_r-(wing.c_r-wing.c_t)*abs(linspace(-1,1,length(xw1)));

zw1 = 0.85*fusRadius*ones(size(xw1)); % height of center of wing
zw2 = wing.t_r-(wing.t_r-wing.t_t)*abs(linspace(-1,1,length(xw1))); % thickness of wing

h(4) = surface(x+[1*xw1;1*xw1;1*xw1],y+[yw1;yw2;yw3],...
    z+[zw1;zw1+.5*zw2;zw1]*zscale,'facecolor',wingColor,...
    'linestyle',linestyle,'edgecolor',edgeColor);
h(5) = surface(x+[1*xw1;xw1;1*xw1],y+[yw1;yw2;yw3],...
    z+[zw1;zw1-.5*zw2;zw1]*zscale,'facecolor',wingColor,...
    'linestyle',linestyle,'edgecolor',edgeColor);

% tail wing (horizontal tail): 
xtw1 = -linspace(-htail.b/2,htail.b/2,10); 
xtw  = [xtw1;xtw1;xtw1];
ytw1 = htail.h_l_r+(htail.h_l_t-htail.h_l_r)*abs(linspace(-1,1,length(xtw1)));  % leading edge of horizontal tail
ytw2 = htail.h_l_r+wing.c_r*0.3+zeros(size(xtw1)); % horitail max thickness line
ytw3 = ytw1+htail.c_r-(htail.c_r-htail.c_t)*abs(linspace(-1,1,length(xtw1))); % horizontal tail trailing edge
ytw = [ytw1;ytw2;ytw3]; 
ztw1 = .9*fusRadius*ones(size(xtw1)); 
h(6) = surface(x+xtw,y+ytw,z+[ztw1;ztw1+.05*fusRadius;ztw1]*zscale,...
    'facecolor',tailWingColor,'linestyle',linestyle,...
    'edgecolor',edgeColor); 
h(7) = surface(x+xtw,y+ytw,z+[ztw1;ztw1;ztw1]*zscale,...
    'facecolor',tailWingColor,'linestyle',linestyle,...
    'edgecolor',edgeColor);

% Fin (vertical tail): 
zts1 = linspace(0,vtail.b,10)+fuse.D/2; 
zts  = [zts1;zts1;zts1];
yts1 = vtail.h_l_r+(vtail.h_l_t-vtail.h_l_r)*abs(linspace(0,1,length(zts1)));  % leading edge of horizontal tail
yts2 = vtail.h_l_r+wing.c_r*0.3+zeros(size(zts1)); % horitail max thickness line
yts3 = yts1+vtail.c_r-(vtail.c_r-vtail.c_t)*abs(linspace(0,1,length(zts1))); % horizontal tail trailing edge
yts  = [yts1; yts2; yts3];
xts1 = vtail.t_r-(vtail.t_r-vtail.t_t)*abs(linspace(-1,1,length(zts1))); % thickness of wing
xts  = [zeros(size(xts1));xts1;zeros(size(xts1))];

h(8) = surface(x+xts,y+yts,z+zts*zscale,...
    'facecolor',finColor,'linestyle',linestyle,...
    'edgecolor',edgeColor);
h(9) = surface(x-xts,y+yts,z+zts*zscale,...
    'facecolor',finColor,'linestyle',linestyle,...
    'edgecolor',edgeColor);

% Propellers: 
xp = prop.D/2*sin(0:.2:2*pi); 
zp = prop.D/2*cos(0:.2:2*pi)+.4*chordlength; %TODO: adjust height
yp = zeros(size(xp))+prop.h; 
h(10) = patch(x+xp,y+yp,z+zp*zscale,propColor); 


% % Propeller Hub: 
% [xph,zph,yph] = cylinder(fusRadius.*([1 .95 .9 .8 .7 .5 .5 .5 .5]).*(cos(linspace(0,pi/2,9)).^.2)); 
% zcn(6:end,:) = zcn(6:end,:)-fusRadius/5; % lower the tip of nose down
% % ycn = -ycn.*.7*chordlength; 
% ycn = -ycn*fusRadius+fuse.D/2;
% h(2) = surface(x+xcn,y+ycn,z+zcn*zscale,...
%     'facecolor',fusColor,'linestyle',linestyle,...
%     'edgecolor',edgeColor);

set(h(10),'facealpha',.2,'edgealpha',.5)


%% Set Roll, Pitch, and Yaw: 

if setroll
    rotate(h,[0 1 0],roll,[x y z]);
end

if setpitch
    rotate(h,[1 0 0],pitch,[x y z]);
end

if setyaw
    rotate(h,[0 0 1],yaw,[x y z]);
end

%% Set view

if SetView
    view([320 30]); 
    xlabel('X (ft)'),ylabel('Y (ft)'),zlabel('Z (ft)');
    axis tight equal
    lighting gouraud
    camlight
else
    axis auto
end

%% Clean up: 

% Return axes to initial hold state
if ~initialHoldState
    hold off
end
    
% Discard surface handles if they're unwanted: 
if nargout==0
    clear h
end
end
