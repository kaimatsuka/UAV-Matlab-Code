%waypoint guidance
%D.Toohey

function  out = my_guid(in)

global traj_total

ind      = in(1);
pE       = in(2);
pN       = in(3);
psi      = in(4);
V        = in(5);
state    = in(6);
% STATE = 1 get to epicenter
%       = 2 spiral 
%       = 3 go and loiter at pt 1
%       = 4 go and loiter at pt 1
%       = 4 go and loiter at pt 1
%       = 99 arrive landing site (catcher)

g = 9.81; %(m/s^2) gravitational acceleration

% Design parameter
L1 = 1000; % ft

% pos=[pE pN]; 
pos=[pN pE];        

% initial iteration
if state == 0
    load_traj;  
    state = 1;
end
    
if state == 1
    traj_remain = traj_total(ind:end,:);
    p_diff = traj_remain-ones(size(traj_remain,1),1)*pos;
    d = (p_diff(:,1).^2 + p_diff(:,2).^2).^0.5;
    ind_remain = find(d>L1);
    
    if isempty(ind_remain)
        ind = size(traj_total,1); % last index of traj total
        state = 99; % move to next point
    else
        ind = ind+ind_remain(1)-1; %update reference location index
    end
end
          

%{
if isempty(ind_remain)
    ind = size(traj_total,1); % last index of traj total
    state = 99; % arrived landing site (catcher)
else
    ind = ind+ind_remain(1)-1; %update reference location index
end
%}

ref_pos = traj_total(ind,:); % reference position

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
out(5) = state;
