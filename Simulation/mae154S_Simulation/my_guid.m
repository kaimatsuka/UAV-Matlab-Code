%waypoint guidance
%D.Toohey

function  out = my_guid(in)

global traj_total

ind      = in(1);
pE       = in(2);
pN       = in(3);
psi      = in(4);
V        = in(5);
V_comm   = in(6);
alt_comm = in(7);
state    = in(8);
% STATE = 1 get to epicenter
%       = 2 spiral 
%       = 3 go and loiter at pt 1
%       = 4 go and loiter at pt 1
%       = 4 go and loiter at pt 1
%       = 99 arrive landing site (catcher)

% constants
g = 9.81; %(m/s^2) gravitational acceleration
mph2fps = 1.4667;

% V_stall  = 40*mph2fps;
% V_loiter = 50*mph2fps;
V_loiter = 100;
V_cruise = 90*mph2fps;

% Lattitude Longitutde of the target
latlon_SFO = [37.621528  -122.376914]; % SFO
latlon_TP  = [37.751937, -122.447410]; % Twin Peaks
latlon_CT  = [37.794242, -122.407179]; % China town
latlon_BH  = [37.734960, -122.412883]; % Bernal Heights
latlon_GGP = [37.768702, -122.462160]; % Golden Gate Park


% lat lon of key locations
latlon_launch = latlon_SFO;
latlon_epic   = latlon_TP;
latlon_loit1  = latlon_CT;
latlon_loit2  = latlon_BH;
latlon_home   = latlon_SFO;
  

% convert lat lon to NED
pos_epic  = latlon2NE(latlon_launch, latlon_epic, latlon_BH(1));
pos_loit1 = latlon2NE(latlon_launch, latlon_loit1, latlon_BH(1));
pos_loit2 = latlon2NE(latlon_launch, latlon_loit2, latlon_BH(1));
pos_home  = [0     0];


% Guidance Gain Parameter
L1 = 1000; % ft

if state == 0
    % initialize trajectory
    traj_total = [0 0];
end
   
% Calculate next point
traj_remain = traj_total(ind:end,:);
pos=[pN pE]; 
p_diff = traj_remain-ones(size(traj_remain,1),1)*pos;
d = (p_diff(:,1).^2 + p_diff(:,2).^2).^0.5;
ind_remain = find(d>L1);

% if we consumed all points on this part of trajectory, load next part of
% trajectory

while isempty(ind_remain)        
    
    ind = 1; % index reset to 1
    state = state+1; % next mode (spiral)
    switch state
        case 1
            traj_total = draw_line(traj_total(end,1),traj_total(end,2),pos_epic(1),pos_epic(2));
            V_comm = V_cruise;
            alt_comm = 7500;
        case 2
            traj_total = draw_spiral(traj_total(end,1),traj_total(end,2));
            V_comm = V_cruise;
            alt_comm = 7500;
        case 3
            traj_total = draw_loiter(traj_total(end,1),traj_total(end,2),pos_loit1(1),pos_loit1(2),pos_loit2(1),pos_loit2(2),3);
            V_comm = V_loiter;
            alt_comm = 1000;
        case 4
            traj_total = draw_loiter(traj_total(end,1),traj_total(end,2),pos_loit2(1),pos_loit2(2),pos_home(1),pos_home(2),3);
            V_comm = V_loiter;
            alt_comm = 1000;
        case 5
            traj_total = draw_line(traj_total(end,1),traj_total(end,2),pos_home(1),pos_home(2));
            V_comm = V_loiter;
            alt_comm = 1000;
    end

    % calculate next point in new trajectory
    traj_remain = traj_total(ind:end,:);
    p_diff = traj_remain-ones(size(traj_remain,1),1)*pos;
    d = (p_diff(:,1).^2 + p_diff(:,2).^2).^0.5;
    ind_remain = find(d>L1);
end
          
if state == 6 %last step
    ref_pos = pos_home;
else
    ind = ind+ind_remain(1)-1; %update reference location index
    ref_pos = traj_total(ind,:); % reference position
end

diff_pos = ref_pos-pos;
ref_head = atan2(diff_pos(2),diff_pos(1)); %(rad) reference location direction
eta = ref_head-psi; %(rad) difference btw heading and reference direction

a_cmd = 2*V^2/L1*sin(eta); %lateral acceleration command

phi_comm = atan(a_cmd/g); 

% banking threashold
phi_max = 30; %(deg)
if phi_comm > phi_max*pi/180
    phi_comm = phi_max*pi/180;
elseif phi_comm < -phi_max*pi/180;
    phi_comm = -phi_max*pi/180;
end

out(1) = ind;
out(2) = phi_comm;
out(3) = ref_pos(2); % reference E
out(4) = ref_pos(1); % reference N
out(5) = V_comm;
out(6) = alt_comm;
out(7) = state;
