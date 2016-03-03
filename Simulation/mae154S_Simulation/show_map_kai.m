%map display function -- D.Toohey
function  out = show_map(in)

pn    = in(1);
pe    = in(2);
alt   = in(3);
tar_E = in(4);
tar_N = in(5);
t     = in(6);

plot(pe,pn,'.b'), hold on

if t==0
    %conversion
    mile2ft = 5280;

    % Image of SF
    img = imread('SF2.jpg');
    % lat lon range of the image
    img_lat = [37.581951 37.831787];
    img_lon = [-122.541591 -122.223880];

    % lat lon of launch location
    launch_lat = 37.621528;
    launch_lon = -122.376914;

    % Earth Radius
    E_R = 3959*mile2ft; 

    % Move [0,0] coordinate of plot to launch location
    y_range = E_R*(img_lat-launch_lat)*pi/180;
%     y_range_inv = [y_range(2) y_range(1)]; % ad hoc solution to plot reversal issue
    y_range = fliplr(y_range); % ad hoc solution to plot reversal issue
    x_range = E_R*cosd(sum(img_lat)/2)*(img_lon-launch_lon)*pi/180;
    figure(1)
    imagesc(x_range, y_range, img);
    set ( gca, 'ydir', 'normal' )
%     error('hold on!')
    hold on
end

plot(tar_E,tar_N,'xr')
xlabel('pE')
ylabel('pN')
axis equal
grid on

%{

figure(1)
subplot(211)

plot(pe,pn,'.b')
hold on
plot(tar_E,tar_N,'xr')
xlabel('pE')
ylabel('pN')
ylim([-500   40000])
xlim([-40000 500])
axis equal
grid on
subplot(212)
hold on
plot(pn,alt,'.r')
xlabel('pN (ft)')
ylabel('Alt (ft)')
xlim([-15000 15000])
ylim([7000 8500])
grid on
%}

out = pn;