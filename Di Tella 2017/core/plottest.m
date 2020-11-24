
figure(1); clf(1);
hold on
vv_x   = v_x/norm(v_x);
vv_xp  = v_xp/norm(v_xp);
vv_xpp = v_xpp/norm(v_xpp);

vv_ah = v_ah/norm(v_ah);
vv_n  = v_n/norm(v_n);
vv_np = v_np/norm(v_np);
vv_ap = Pfun(vv_n,zeros(3,1),vv_ah);

tmp = [(vv_ap(1)-vv_ah(1))
       (vv_ap(2)-vv_ah(2))
       (vv_ap(3)-vv_ah(3))];

mt = norm(tmp);
t = linspace(0,mt,100);
vx = vv_ah(1)+t*vv_n(1);
vy = vv_ah(2)+t*vv_n(2);
vz = vv_ah(3)+t*vv_n(3);

lh = plot3(vx,vy,vz,'r','Linewidth',1);
% ph = plot3(vx(end),vy(end),vz(end),'*r');

xx = linspace(-2,2,2);
yy = linspace(-2,2,2);
zz = linspace(-2,2,2);

if vv_n(3)~=0
    [x,y] = meshgrid(xx,yy);
    z = (-x*vv_n(1) -y*vv_n(2))/vv_n(3);
else
    [x,z] = meshgrid(xx,zz);
    y = -x*vv_n(1)/vv_n(2);
end
surf(x,y,z)

if vv_np(3)~=0
    [x,y] = meshgrid(xx,yy);
    z = (-x*vv_np(1) -y*vv_np(2))/vv_np(3);
else
    [x,z] = meshgrid(xx,zz);
    y = -x*vv_np(1)/vv_np(2);
end
surf(x,y,z)

v_npp = cross(v_xp,v_xpp);
vv_npp = v_npp/norm(v_np);
if vv_npp(3)~=0
    [x,y] = meshgrid(xx,yy);
    z = (-x*vv_npp(1) -y*vv_npp(2))/vv_npp(3);
else
    [x,z] = meshgrid(xx,zz);
    y = -x*vv_npp(1)/vv_npp(2);
end
surf(x,y,z)

plot3(vv_ap(1),vv_ap(2),vv_ap(3),'or')
plot3(vv_ah(1),vv_ah(2),vv_ah(3),'*b')

plot3(vv_x(1),vv_x(2),vv_x(3),'ok')
plot3(vv_xp(1),vv_xp(2),vv_xp(3),'ok')
plot3(vv_xpp(1),vv_xpp(2),vv_xpp(3),'ok')

Cone([0 0 0],[1/2 0 1/2],[0 sqrt(1/2)],100,'b',1,0);
alpha(0.1)

view([120 40])
axis equal
grid on

xlim([0 1.25])
ylim([-1 1])
zlim([min([0 vv_ap(3)]) 1.25])

xlabel('x')
ylabel('y')
zlabel('z')

pause