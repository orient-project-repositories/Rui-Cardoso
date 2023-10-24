function [lin_ss] = get_linearss(v_val, delta_t, alpha1, alpha2)
%lin_ss System linearization function. Valid only for small theta values.

v=sym('v',[1 3],'real');
w=sym('w',[1 3],'real');
u=sym('u',[1 3],'real');
t=sym('t',1);


%theta ->0
f = [(1/2)*(12-norm(v)^2)/6.*w+cross(w,v)-dot(w,v)*(60+norm(v)^2)/360.*v];
f11=jacobian(f,[v, w]);



%% df2/dv = d(J^-1*tau_ela)/dv
Q1_0 = [-0.007; 0.05302; 0.05302];
Q2_0 = [-0.007; 0.05302; -0.05302];

Q3_0 = [-0.007; -0.07498; 0];
Q4_0 = [-0.007; 0.07498; 0];

Q5_0 = [-0.007; -0.05302; 0.05302];
Q6_0 = [-0.007; -0.05302; -0.05302];

Q_0=[Q1_0 Q2_0 Q3_0 Q4_0 Q5_0 Q6_0];

k = 6;

l1_0 = 34 * 10^-2;
l2_0 = 36 * 10^-2;
l3_0 = 24 * 10^-2;
l4_0 = 24 * 10^-2;
l5_0 = 34 * 10^-2;
l6_0 = 36 * 10^-2;

l_0 = [l1_0 l2_0 l3_0 l4_0 l5_0 l6_0];

% get new P
P1 = [-0.436; -0.0835; 0.021 + 0.1];
P2 = [-0.436; -0.0835; 0.021 - 0.1];

% P3 = [-0.323; -0.140; 0];%-0.323
% P4 = [-0.323; 0.140; 0];%-0.323
P3 = [-0.323; -0.140; 0.0475];%-0.323
P4 = [-0.323; 0.140; 0.0475];%-0.323

P5 = [-0.436; 0.0835; 0.021 + 0.1];
P6 = [-0.436; 0.0835; 0.021 - 0.1];

% concatenate
P = [P1 P2 P3 P4 P5 P6];

%----------------------------
% connection points
X1 = [-0.2055; 0; 0.0525];
X2 = [-0.2055; 0; -0.0525];
X5 = X1;
X6 = X2;

% intermediate string direction points
Yi = [X1 X2 P(:,3) P(:,4) X5 X6];

Q=sym('Q',[3 6]);
Xi = [X1 X2 Q(:,3) Q(:,4) X5 X6];


tau_ela=sym('tau_ela',[3 6]);
dtau_eladQ=sym('dtau_eladQ', [3 3 3]);

dtau_eladv=sym(zeros(3));
df2dv=sym(zeros(3));


inertia=[ 0.0004759 0 0;0 0.0004759 0;0 0 0.0004759];


dRdv=sym('dRdv',[3 3 3]);

for i=1:6
    l=norm(P(:,i)-Xi(:,i))+norm(Xi(:,i)-Q(:,i));
    d=(Yi(:,i)-Q(:,i))/norm(Yi(:,i)-Q(:,i));
    F_v(:,i)=(l-l_0(i))*d;
    F_v(:,i)=k/l_0(i)*F_v(:,i);
    tau_ela(:,i)=cross(Q(:,i),F_v(:,i));
    dtau_eladQ(:,:,i)=jacobian(tau_ela(:,i),Q(:,i));
    
    for j=1:3
        dtau_eladv(:,j)=dtau_eladv(:,j)+dtau_eladQ(:,:,i)*dRdv(:,:,j)*Q_0(:,i);
        df2dv(:,j)=df2dv(:,j)+inertia\(dtau_eladQ(:,:,i)*dRdv(:,:,j)*Q_0(:,i));
    end
   
    
    
end



R=expm([0 -v_val(3) v_val(2); v_val(3) 0 -v_val(1); -v_val(2) v_val(1) 0]);

dtau_gdv=get_dtau_gdv(v_val);

for i=1:3
    dtau_eladv(:,i)=dtau_eladv(:,i)+dtau_gdv(:,i);
end


df2dv=subs(df2dv, Q, R*Q_0); 
df2dv=subs(df2dv,dRdv,get_dRdv(v_val));

df2dw=jacobian(-inertia^-1*0.02*w',w);
df2du=jacobian(inv(inertia)*u'/alpha2,u);

A = [f11; df2dv df2dw];
B = [zeros(3);df2du];

B = subs(B, t, delta_t*t);


A=subs(A,[v w],[v_val v_val]);
H=[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0];

lin_ss=ss(double(A),double(B),H,0);
lin_ss=c2d(lin_ss, delta_t);
end

