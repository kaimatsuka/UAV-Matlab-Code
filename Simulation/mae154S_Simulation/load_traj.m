function [] = load_traj()

global traj_total

    % Lattitude Longitutde of the target
    
    latlon_SFO = [37.621528  -122.376914]; % SFO
    latlon_TP  = [37.751937, -122.447410]; % Twin Peaks
    latlon_CT  = [37.794242, -122.407179]; % China town
    latlon_BH  = [37.734960, -122.412883]; % Bernal Heights
    latlon_GGP = [37.768702, -122.462160]; % Golden Gate Park

    
    % lat lon
    latlon_launch = latlon_SFO;
    latlon_epic   = latlon_TP;
    latlon_loit1  = latlon_CT;
    latlon_loit2  = latlon_BH;
    latlon_home   = latlon_SFO;
  
    
    % North East location of target positinos relative to launch position
%     pos_epic  = [50000 -17000];
%     pos_loit1 = [55000 -10000];
%     pos_loit2 = [55000 -20000];
%     pos_home  = [0     0];
    
    pos_epic  = latlon2NE(latlon_home, latlon_epic, latlon_BH(1));
    pos_loit1 = latlon2NE(latlon_home, latlon_loit1, latlon_BH(1));
    pos_loit2 = latlon2NE(latlon_home, latlon_loit2, latlon_BH(1));
    pos_home  = [0     0];
    
    
    [traj1 heading] = draw_line(0,-300,pos_epic(1),pos_epic(2));
    [traj2 heading] = draw_spiral(traj1(end,1),traj1(end,2));
    [traj3 heading] = draw_loiter(traj2(end,1),traj2(end,2),pos_loit1(1),...
                                 pos_loit1(2),pos_loit2(1),pos_loit2(2),3);
    [traj4 heading] = draw_loiter(traj3(end,1),traj3(end,2),pos_loit2(1),...
                                 pos_loit2(2),pos_home(1),pos_home(2),3);
    [traj5 heading] =  draw_line(traj4(end,1),traj4(end,2),pos_home(1),pos_home(2));
    
    traj_total = traj1;
%     traj_total = [traj1;
%                   traj2;
%                   traj3;
%                   traj4;
%                   traj5];
%         
end      

