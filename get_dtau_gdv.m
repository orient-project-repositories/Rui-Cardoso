function [dtau_gdv] = get_dtaug_dv(v)
%Partial derivative of gravity torque w.r.t angular position v

dtau_gdv=sym(zeros(3));
R=expm([0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]);
dRdv=get_dRdv(v);
center_mass = R*[0.02; 0; 0]; % center of mass
c_mass=sym('c_mass',[3 1]);
mass = 0.248; % mass in kg
Fg = mass*[0; 0; -9.8];
tau_g = cross(c_mass, Fg);

dtau_gdc_mass=jacobian(tau_g,c_mass);

for i=1:3
    dtau_gdv(:,i)=dtau_gdc_mass*dRdv(:,:,i)*center_mass;
end

end

