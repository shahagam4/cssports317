function dy=projectile_fn(t,y)
    dy=zeros(7,1);
    global mass Cd air_density Wx Wy Wz W Cl Area;
    dy(1)=y(4);
    dy(2)=y(5);
    dy(3)=y(6);
    dy(4)=-(0.5*air_density*Area*sqrt(y(4)^2+y(5)^2+y(6)^2)*(Cd*y(4)-Cl*((Wy*y(6)-Wz*y(5))/W)))/mass;
    dy(5)=-(0.5*air_density*Area*sqrt(y(4)^2+y(5)^2+y(6)^2)*(Cd*y(5)-Cl*((Wz*y(4)-Wx*y(6))/W)))/mass;
    dy(6)=-9.8-(0.5*air_density*Area*sqrt(y(4)^2+y(5)^2+y(6)^2)*(Cd*y(6)-Cl*((Wx*y(5)-Wy*y(4))/W)))/mass;
    dy(7) = 0.5*air_density*Area*Cd*sqrt(y(4)^2+y(5)^2+y(6)^2)^3;
end