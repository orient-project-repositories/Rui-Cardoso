%% SDN
% 
% 
function [G,K,opt_p,mincost] = todorov_getgains_p(ssmodelo, qx, qy, C, D, s, x_init,L)
%todorov_getgains - computation of controller and kalman gains using emo
%todorov's iterative method
%   Detailed explanation goes here
%set_param('ssmodelo_kf_getgains','FastRestart','off');

 lambda_endp=1;%.1;%1e10
 lambda_est_e=1;%1e10
 lambda_t=1;%1e-3%e-4;%0;

alpha_reward=150e3;%.4e18
beta_reward=.095;%.01;%3

Q_x = qx*eye(size(ssmodelo.A,1));
Q_y = qy*eye(size(ssmodelo.C,1));

k1=20;

A=ssmodelo.A; 
B=ssmodelo.B;
H=ssmodelo.C;
%H=eye(6);
delta_t=ssmodelo.Ts;
clear ssmodelo;

%D=10*D;
%C=50*C;

L_shadmehr = 0*.5*[1 0 0; 0 1 0; 0 0 1];
T=150*[0 0 0 0 0 0; 0 1e6 0 0 0 0; 0 0 1e6 0 0 0; 0 0 0 1.5e1 0 0; 0 0 0 0 1e2 0; 0 0 0 0 0 1e2]; %AD linearized sys
T_during=[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 .1 0; 0 0 0 0 0 .1]; 

%T=[5e9 0 0 0 0 0; 0 5e9 0 0 0 0; 0 0 5e9 0 0 0; 0 0 0 1e6 0 0; 0 0 0 0 1e6 0; 0 0 0 0 0 1e6]; %shadmehr values
%T=H'*[10e3 0 0; 0 10e3 0; 0 0 10e3]*H;
%T=[.01*1e9 0 0 0 0 0; 0 5e9 0 0 0 0; 0 0 5e9 0 0 0; 0 0 0 0.01*1e6 0 0; 0 0 0 0 1e5 0; 0 0 0 0 0 1e5]; 

%T=[10 0 0 0 0 0; 0 10 0 0 0 0; 0 0 10 0 0 0; 0 0 0 .00001 0 0; 0 0 0 0 .00001 0; 0 0 0 0 0 .00001]; %LP-AD
% T=[10e20 0 0; 0 10e20 0; 0 0 10e20]; %identified sys
%T = [10 0 0; 0 0.00001 0; 0 0 0.1]; %1Dmodel

n = size(A,1);
m = size(B,2);
out_size=size(H,1);

mincost=inf;
pmax=200;
pmin=5;

pmin=30
pmax=30
step_p_search=2;

for p=pmin:step_p_search:pmax

     L_shadmehr=L*(.0001/(p*delta_t))*eye(m);
     
G = zeros(m,n,p); %controller gain
P = zeros(n,n,p);   %prediction kalman
K = zeros(n,out_size,p);   %kalman gain

 
W_x = zeros(n,n,p);
W_e = zeros(n,n,p);
w = zeros(1,p);
 
C_x = zeros(m,m,p);
C_e = zeros(m,m,p);
D_e = zeros(n,n,p);

 %starting from the last point - W_x = T, W_e = 0
W_x(:,:,p)=T;
W_e(:,:,p)=zeros(n);
w(p)=(alpha_reward*beta_reward/(1+beta_reward*p));
P(:,:,1)=[pi/4 0 0 0 0 0; 0 (pi/4) 0 0 0 0; 0 0 pi/4 0 0 0; 0 0 0 (pi/4)*3 0 0; 0 0 0 0 (pi/4)*3 0; 0 0 0 0 0 (pi/4)*3];

 %one begins by computing a set of Kalman gains K .These gains can be computed 
    %from Eqs. (12.37) and (12.38) with the signal dependent noise component set to zero
for k=1:p
    

    K(:,:,k)=P(:,:,k)*H'*(H*P(:,:,k)*H'+Q_y)^-1; % dimensao de Q_y
    P(:,:,k)=P(:,:,k)*(eye(n)-H'*K(:,:,k)');
    P(:,:,k+1)=A*P(:,:,k)*A'+Q_x;
end



for it=1:3
    %One then computes a set of control gains G using Eqs. (12.33) and (12.35)
    for k=p-1:-1:1
        for i=1:m
            C_x(:,:,k+1)=C_x(:,:,k+1)+C(:,:,i)'*B'*W_x(:,:,k+1)*B*C(:,:,i);
            C_e(:,:,k+1)=C_e(:,:,k+1)+C(:,:,i)'*B'*W_e(:,:,k+1)*B*C(:,:,i);
        end
        for i=1:n
            D_e(:,:,k+1)=D_e(:,:,k+1)+D(:,:,i)'*H'*K(:,:,k)'*A'*W_e(:,:,k+1)*A*K(:,:,k)*H*D(:,:,i);
        end


    %nao é preciso calcular wx we nem w na 1a iteração
    % G nao depende do reward? (-w(k+1))
        G(:,:,k)=((L_shadmehr+C_x(:,:,k+1)+C_e(:,:,k+1)+B'*W_x(:,:,k+1)*B)^-1)*B'*W_x(:,:,k+1)*A;

        W_e(:,:,k)=((A-A*K(:,:,k)*H)'*W_e(:,:,k+1)*(A-A*K(:,:,k)*H)+G(:,:,k)'*B'*W_x(:,:,k+1)*A);%lambda_est_e*
        W_x(:,:,k)=(A'*W_x(:,:,k+1)*A+D_e(:,:,k+1)-G(:,:,k)'*B'*W_x(:,:,k+1)*A);%lambda_endp*
        w(k)= (trace(W_x(:,:,k+1)*Q_x+W_e(:,:,k+1)*(Q_x+A*K(:,:,k)*Q_y*K(:,:,k)'*A')))+w(p);%lambda_t*

    end


    %One then re-computes the Kalman gains using Eq. (12.39)
    S_e(:,:,1)=P(:,:,1); %??
    S_x(:,:,1)=eye(n);
    for k=1:p
         %K(:,:,k)=eye(6);
        
        %shadmehr errado-----------
        D_aux=0;
        for i=1:n
            D_aux = D_aux + H*D(:,:,i)*(S_x(:,:,k))*D(:,:,i)'*H';
        end
        %--------------------------
        K(:,:,k)=S_e(:,:,k)*H'*(H*S_e(:,:,k)*H'+Q_y+D_aux)^-1;
        if k<p
            C_aux=0;
            for i=1:m
                C_aux = C_aux + B*C(:,:,i)*G(:,:,k)*S_x(:,:,k)*G(:,:,k)'*C(:,:,i)'*B';
            end
            S_e(:,:,k+1)=Q_x+(A-A*K(:,:,k)*H)*S_e(:,:,k)+C_aux;
            S_x(:,:,k+1)=A*K(:,:,k)*H*S_e(:,:,k)*A'+(A-B*G(:,:,k))*S_x(:,:,k)*(A-B*G(:,:,k))';
        end
    end
    
    
end

 %simular para calcular o custo total
 %options = simset('SrcWorkspace','current');
 %simout=sim('ssmodelo_kf_getgains',[],options);
 %sim('olho_kf_getgains',[],options); %1D
 index=(p-pmin)/step_p_search+1;
 cost(index)=0;
 cost_t(index)=0;
 cost_endp_a(index)=0;
 cost_est_e(index)=0;
 

 
 %est=simout.est;
 x=zeros(n,p);
 u=zeros(m,p);
 y=zeros(out_size,p);
 est=zeros(n,p);
 s_x=zeros(n,1);
 

 
 for k=1:p+1
     %calculo do custo para um dado p
     %xk=simout.state_x(k,:)';
        
     if k==1
         u=zeros(size(B,2),1);
         x_prev=x_init;
         est=zeros(size(A,1),1);
     else
         u=G(:,:,k-1)*(s-x_prev);
%     s_x=A*s_x+B*uff;
        x(:,k-1)=A*x_prev+B*(u);%+normrnd(0,qx^2,n,1);
        y(:,k-1)=H*x(:,k-1);%+normrnd(0,qy^2,m,1);
     
        est=A*est+A*K(:,:,k-1)*(y(:,k-1)-H*est)+B*u;
        x_prev=x(:,k-1);
     
     
     
     
     
   
    xk=x(:,k-1);
    cost(index)=cost(index)+lambda_endp*abs(s-xk)'*W_x(:,:,k-1)*abs(s-xk)+lambda_est_e*(xk-est)'*W_e(:,:,k-1)*(xk-est)+lambda_t*w(k-1); 
     cost_t(index)=cost_t(index)+lambda_t*w(k-1);
     cost_endp_a(index)=cost_endp_a(index)+lambda_endp*abs(s-xk)'*W_x(:,:,k-1)*abs(s-xk);
     cost_est_e(index)=cost_est_e(index)+lambda_est_e*(xk-est)'*W_e(:,:,k-1)*(xk-est);
%     cost(index)=cost(index)+abs(s-xk)'*W_x(:,:,k-1)*abs(s-xk)+(xk-est)'*W_e(:,:,k-1)*(xk-est)+w(k-1); 
%      cost_t(index)=cost_t(index)+w(k-1);
%      cost_endp_a(index)=cost_endp_a(index)+abs(s-xk)'*W_x(:,:,k-1)*abs(s-xk);
%      cost_est_e(index)=cost_est_e(index)+(xk-est)'*W_e(:,:,k-1)*(xk-est);
     end
 end
    cost;

  %  bla
 if cost(index)<mincost
     mincost=cost(index);
     opt_p=p;
     G_out=G;
     K_out=K;
 end
         
 
end

G=G_out;
K=K_out;

p_array=pmin:step_p_search:pmax;

end

