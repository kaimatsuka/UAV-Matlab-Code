%waypoint guidance
%D.Toohey

function  out = wayguid(in)

way_num = in(1);
pE = in(2);
pN = in(3);
heading = in(4);

if way_num == 0
    way_num = 1;
end

%pre-defined waypoints  (pE   pN
waypoints = [0      36000;...
             -3000  36000;...
             -3000  0;...
             -6000  0;...
             -6000  36000;...
             -9000  36000;...
             -9000  0;...
             -12000 0;...
             -12000 36000;...
             -15000 36000;...
             -15000 0;...
             -18000 0;...
             -18000 36000;...
             -21000 36000;...
             -21000 0;...
             -24000 0;...
             -24000 36000;...
             -27000 36000];
         
% distance threshold used to switch to next waypoint
dist_thresh = 50;
         
tar_E = waypoints(way_num,1);
tar_N = waypoints(way_num,2);

delta_E = tar_E - pE;
delta_N = tar_N - pN;
way_dist = (delta_E^2 + delta_N^2)^.5;
if way_dist < dist_thresh
    way_num = way_num + 1;
end

tar_head = atan2(delta_E,delta_N);

delta_psi = tar_head - heading;


%check for angles larger than 180 deg
if delta_psi > pi
    delta_psi = delta_psi - 2*pi;
elseif delta_psi < -pi
    delta_psi = delta_psi + 2*pi;
end

phi_comm = .6*delta_psi;

if phi_comm > 30*pi/180
    phi_comm = 30*pi/180;
elseif phi_comm < -30*pi/180;
    phi_comm = -30*pi/180;
end

out(1) = way_num;
out(2) = phi_comm;
out(3) = tar_E;
out(4) = tar_N;
         
         
   
   
       
