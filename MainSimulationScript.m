clear;
clc;
digits(100)
format long
warning('off', 'MATLAB:nearlySingularMatrix')
nr=15; %number of nodes in radial direction --> Main mesh 
nz=15; % number of nodes in z direction 
dt=0.01;  % time step in HC equation 
a=100e-3;  % inner Radious of plate
b=1e-8;  % outer Radious of plate
h=5e-3; % thickness m  h is thickness of circular plate
kisi=0;
k_shear=5/6;
% 'insul' --> zero heat flux 
% 'convection' --> heat convection with air
% 'prescribed' --> prescribed temprature
right_edge='insul';
top_edge='prescribed';
bottom_edge='prescribed';
left_edge='insul';

partial_ratio_top=0.55; % percetatge of the partial load on top face
partial_ratio_bottom=1; % percetatge of the partial load on bottom face

[c0,c1,c2,r]= GDQ_1D(nr,a,b);
[TT,tt,nt,r,zz,kappa,cp,ro, E, alpha_thermal, nu,index_heated_top,all_index_top]=general_heat_two_dimension_FGM_GDQ_TID(nr,nz,dt,a,b,h,kisi,right_edge,top_edge,bottom_edge,left_edge,partial_ratio_top,partial_ratio_bottom);

%%
% Numerical  Trapezoidal  integration for A, B, D matrix, N_T & M_T, I1, I2, I3 

k=zeros(3*nr,3*nr,nt);
F=zeros(3*nr,1,nt);

    % A B D matrixes would be "r" dependent in case of 2D FGM
 for i=1:nt
     for j=1:nr
    N_TT(j,1,i) = trapz(-zz,E(:,j).*alpha_thermal(:,j).*(TT(:,j,i)-TT(:,j,1))./(1-nu(:,j))); % numerical integration with trapeodial method   % N_TT is dependent on the "r" in cas eof axussymetric partial load 
    M_TT(j,1,i) = trapz(-zz,zz'.*E(:,j).*alpha_thermal(:,j).*(TT(:,j,i)-TT(:,j,1))./(1-nu(:,j))); % numerical integration with trapeodial method  % M_TT is dependent on the "r" in cas eof axussymetric partial load 
    
    A11(j,1,i)=trapz(-zz,E(:,j)./(1-nu(:,j).^2));
    A22(j,1,i)=trapz(-zz,E(:,j)./(1-nu(:,j).^2));
%A22=trapz(-zz,E./(1-nu.^2));
    A12(j,1,i)=trapz(-zz,(nu(:,j).*E(:,j))./(1-nu(:,j).^2));
    A55(j,1,i)=(k_shear)*trapz(-zz,E(:,j)./(2.*(1+nu(:,j))));
    B11(j,1,i)=trapz(-zz,(zz'.*E(:,j))./(1-nu(:,j).^2));
    B22(j,1,i)=trapz(-zz,(zz'.*E(:,j))./(1-nu(:,j).^2));
%B22=trapz(-zz,(zz'.*E)./(1-nu.^2));
    B12(j,1,i)=trapz(-zz,(zz'.*nu(:,j).*E(:,j))./(1-nu(:,j).^2));
    D11(j,1,i)=trapz(-zz,(((zz').^2).*E(:,j))./(1-nu(:,j).^2));
    D12(j,1,i)=trapz(-zz,(((zz').^2).*nu(:,j).*E(:,j))./(1-nu(:,j).^2));
    D22(j,1,i)=trapz(-zz,(((zz').^2).*E(:,j))./(1-nu(:,j).^2));

     end
 end
 

% figure(1)
% hold on 
% plot(r,reshape(N_TT(:,1,end),[length(r),1]));
% 
% figure(2)
% hold on 
% plot(r,reshape(M_TT(:,1,end),[length(r),1]));
    
    %%
% Linear Stiffness
for i=1:nt
    
k(1:nr,1:nr,i)=(A11(:,1,i)).*(c2+(r.^-1)'.*c1 -((r.^2).^-1)'.*c0); %K_uu
k(1:nr,nr+1:2*nr,i)=0; %K_uw
k(1:nr,2*nr+1:3*nr,i)=(B11(:,1,i)).*(c2 +(r.^-1)'.*c1 -((r.^2).^-1)'.*c0); %K_u phi

k(nr+1:2*nr,1:nr,i)=0; %K_wu
k(nr+1:2*nr,nr+1:2*nr,i)=(A55(:,1,i)).*(c2 +(r.^-1)'.*c1) -(c0*N_TT(:,1,i)).*(c2 +(r.^-1)'.*c1)  -(c1*N_TT(:,1,i)).*(c1);   %K_ww  % N_TT is changin in "r" direction 
k(nr+1:2*nr,2*nr+1:3*nr,i)=(A55(:,1,i)).*(c1 +(r.^-1)'.*c0);%K_w phi

k(2*nr+1:3*nr,1:nr,i)=(B11(:,1,i)).*(c2 +(r.^-1)'.*c1 -((r.^2).^-1)'.*c0); % K_phi u 
k(2*nr+1:3*nr,nr+1:2*nr,i)=-A55(:,1,i).*c1;% % K_phi W 
k(2*nr+1:3*nr,2*nr+1:3*nr,i)=(D11(:,1,i)).*(c2 +(r.^-1)'.*c1 -((r.^2).^-1)'.*c0) -A55(:,1,i).*(c0);% % K_phi Phi 

%BC
% always at r=0 --> u=0 & Q_rz + N_rr*w,r=0 & phi=0  
% Full circular IM S

for J=1:nr
 
    F(J,1,i)=c1(J,:)*N_TT(:,1,i);   % effect of partial load d(N_TT)/dr on the RHS 

end

for JJ=2*nr+1:3*nr
 
    F(JJ,1,i)=c1(JJ-2*nr,:)*M_TT(:,1,i);   % effect of partial load d(M_TT)/dr on the RHS 

end


k(1,:,i)=0; 
k(1,1,i)=1; 
F(1,1,i)=0;

k(nr,:,i)=0;
k(nr,nr,i)=1;
F(nr,1,i)=0;


k(nr+1,1:nr,i)=0; % for r=0 and Qrz +Nrr*w,r ==0 (linear terms are all zero becuase of r=0 multiplier ) 
k(nr+1,nr+1:2*nr,i)=A55(1,1,i).*c1(1,:) -(c0(1,:)*N_TT(:,1,i))*c1(1,:);
k(nr+1,2*nr+1:3*nr,i)=A55(1,1,i).*c0(1,:);

k(2*nr,:,i)=0;
k(2*nr,2*nr,i)=1;

k(2*nr+1,:,i)=0;
k(2*nr+1,2*nr+1,i)=1;
F(2*nr+1,1,i)=0;

k(3*nr,1:nr,i)=B11(nr,1,i).*c1(nr,:) +B12(nr,1,i).*c0(nr,:)/a;  % r=a in here since Mrr @ r=a is zero
k(3*nr,nr+1:2*nr,i)=0;    % r=a in here since Mrr @ r=a is zero
k(3*nr,2*nr+1:3*nr,i)=D11(nr,1,i).*c1(nr,:) +D12(nr,1,i)*c0(nr,:)/a;   % r=a in here since Mrr @ r=a is zero
F(3*nr,1,i)=c0(nr,:)*M_TT(:,1,i);




UU(:,1,i)=inv(k(:,:,i))*F(:,1,i); % linear solution of [K]*{X}={F}
%W(1,i)=UU((3*nr+1)/2,1,i)/h;
W_linear(:,1,i)=UU(nr+1:2*nr,1,i);
W1_linear(i,1)=W_linear(1,1,i);
Wmid_lin(i,1)=W_linear(ceil(nr/2),1,i);
end

% figure(1)
% %plot(tt,W1_linear/h)
% hold on
% plot(tt,W1_linear/h)
% 
% figure(2)
% hold on 
% plot(r,reshape(W_linear(:,1,end),[nr,1,1])/h)

% start of nonlinear %%  Newton -- Raphson 
u(:,1,1)=zeros(nr,1);
w(:,1,1)=zeros(nr,1);
phi(:,1,1)=zeros(nr,1);
qq(1:3*nr,1,1)=UU(:,1,1);
   
K_linear=k;
for i=2:nt
      
   itera=1;
   sum1=1;
   sum2=1;    

   if i==2 
         
           qq(1:3*nr,1,2)=UU(:,1,2);
           u(:,1,i)=UU(1:nr,1,2);
           w(:,1,i)=UU(nr+1:2*nr,1,2);
           phi(:,1,i)=UU(2*nr+1:3*nr,1,2);
           YY=qq(1:3*nr,1,2);
   else
           u(:,1,i)=u(:,1,i-1);
           w(:,1,i)=w(:,1,i-1);
           phi(:,1,i)=phi(:,1,i-1);
           qq(:,1,i)=[u(:,1,i);w(:,1,i);phi(:,1,i)];
   end    
   
   while sqrt(sum1/sum2)>.0000001
   

        F_linear(:,1,i)=K_linear(:,:,i)*qq(:,1,i) -F(:,1,i); % linear Eq with all applied BC  ###########IMPORTANT "-F" term is for the partial load consideration 
                    
        % nonlinear Equations 
        FNL(1:nr,1,i)=A11(:,1,i).*((c1*w(:,1,i)).*(c2*w(:,1,i)) + ((2*r).^-1)'.*(c1*w(:,1,i)).^2) - A12(:,1,i).*((2*r).^-1)'.*((c1*w(:,1,i)).^2);
        FNL(nr+1:2*nr,1,i)=A11(:,1,i).*(((r).^-1)'.*(c1*u(:,1,i)).*(c1*w(:,1,i)) +((2*r).^-1)'.*(c1*w(:,1,i)).^3  +(c2*w(:,1,i)).*(c1*u(:,1,i)) + (3/2)*(c2*w(:,1,i)).*(c1*w(:,1,i)).^2  +(c1*w(:,1,i)).*(c2*u(:,1,i)))+...
        +A12(:,1,i).*(((r).^-1)'.*(c2*w(:,1,i)).*(c0*u(:,1,i)) + ((r).^-1)'.*(c1*w(:,1,i)).*(c1*u(:,1,i))) +B11(:,1,i).*(((r).^-1)'.*(c1*w(:,1,i)).*(c1*phi(:,1,i)) +(c2*w(:,1,i)).*(c1*phi(:,1,i)) +(c1*w(:,1,i)).*(c2*phi(:,1,i)))+...
        +B12(:,1,i).*(((r.^2).^-1)'.*(c1*w(:,1,i)).*(c0*phi(:,1,i)) +((r).^-1)'.*(c1*w(:,1,i)).*(c1*phi(:,1,i)));
        FNL(2*nr+1:3*nr,1,i)=B11(:,1,i).*((c1*w(:,1,i)).*(c2*w(:,1,i)) +((2*r).^-1)'.*(c1*w(:,1,i)).^2) -B12(:,1,i).*((2*r).^-1)'.*((c1*w(:,1,i)).^2);

          
        FFF(:,1,i)=F_linear(:,1,i)+FNL(:,1,i); 
        
        % kole moadelate system shamele khati va gheire khati ba ehtesabe
        % sharte marzi -->  % boundry condition full circular IM S for Initial terms  (only effects on
        % FNL and Flinear) 
        
        FFF(1,1,i)=0;% for full circular
        FFF(nr,1,i)=0;
        FFF(nr+1,1,i)=A55(1,1,i).*(c1(1,:)*w(:,1,i) +c0(1,:)*phi(:,1,i)) +(A11(1,1,i).*c1(1,:)*u(:,1,i) +A11(1,1,i).*0.5*(c1(1,:)*w(:,1,i)).^2 +A12(1,1,i).*c0(1,:)*u(:,1,i)/b +B11(1,1,i).*c1(1,:)*phi(:,1,i) +B12(1,1,i).*c0(1,:)*phi(:,1,i)/b  -c0(1,:)*N_TT(:,1,i))*(c1(1,:)*w(:,1,i)); % for full circular 
        FFF(2*nr,1,i)=0;
        FFF(2*nr+1,1,i)=0; % for full circular
        FFF(3*nr,1,i)= -c0(nr,:)*M_TT(:,1,i) +B11(nr,1,i).*0.5*(c1(nr,:)*w(:,1,i)).^2 +B11(nr,1,i).*(c1(nr,:)*u(:,1,i)) +B12(nr,1,i).*(c0(nr,:)*u(:,1,i))/a +D11(nr,1,i).*c1(nr,:)*phi(:,1,i) +D12(nr,1,i).*(c0(nr,:)*phi(:,1,i))/a ;
        
        FF_total=FFF(:,1,i);
        
        % differentiotion with respect of each DOF
        
        f11=A11(:,1,i).*(c2 +((r).^-1)'.*(c1) -((r.^2).^-1)'.*(c0));
        f12=A11(:,1,i).*(c1.*(c2*w(:,1,i)) +(c1*w(:,1,i)).*(c2) +((2*r).^-1)'.*((c1).*(c1*w(:,1,i)) +(c1*w(:,1,i)).*(c1))) -A12(:,1,i).*(((2*r).^-1)'.*((c1).*(c1*w(:,1,i))) +((2*r).^-1)'.*((c1*w(:,1,i)).*(c1))); % test shavad
        f13=B11(:,1,i).*(c2 +((r).^-1)'.*(c1) -((r.^2).^-1)'.*(c0));
        
        f21=A11(:,1,i).*(((r).^-1)'.*(c1).*(c1*w(:,1,i)) +(c2*w(:,1,i)).*(c1) +(c1*w(:,1,i)).*(c2)) +A12(:,1,i).*(((r).^-1)'.*(c2*w(:,1,i)).*(c0) +((r).^-1)'.*(c1*w(:,1,i)).*(c1));
        f22=A55(:,1,i).*(c2 +((r).^-1)'.*(c1)) +A11(:,1,i).*(((r).^-1)'.*(c1*u(:,1,i)).*(c1) +3*((2*r).^-1)'.*((c1).^3).*(w(:,1,i)).^2 +(c2).*(c1*u(:,1,i)) +1.5*(c2.*(c1*w(:,1,i)).^2 +(c2*w(:,1,i)).*(2*(c1.^2).*w(:,1,i))) +(c1).*(c2*u(:,1,i)))+...
          +A12(:,1,i).*(((r).^-1)'.*(c2).*(c0*u(:,1,i)) +((r).^-1)'.*(c1).*(c1*u(:,1,i))) +B11(:,1,i).*(((r).^-1)'.*(c1).*(c1*phi(:,1,i)) +(c2).*(c1*phi(:,1,i)) +(c1).*(c2*phi(:,1,i))) +B12(:,1,i).*(((r.^2).^-1)'.*(c1).*(c0*phi(:,1,i)) +((r).^-1)'.*(c1).*(c1*phi(:,1,i))) -(c0*N_TT(:,1,i)).*(c2 +(r.^-1)'.*c1)  -(c1*N_TT(:,1,i)).*(c1);
        f23=A55(:,1,i).*(c1+((r).^-1)'.*c0) +B11(:,1,i).*(((r).^-1)'.*(c1*w(:,1,i)).*(c1) +(c2*w(:,1,i)).*(c1) +(c1*w(:,1,i)).*(c2)) +B12(:,1,i).*(((r.^2).^-1)'.*(c1*w(:,1,i)).*(c0) +((r).^-1)'.*(c1*w(:,1,i)).*(c1));
        
        f31=B11(:,1,i).*(c2 +((r).^-1)'.*c1 -((r.^2).^-1)'.*c0);
        f32=B11(:,1,i).*(c1.*(c2*w(:,1,i)) +(c1*w(:,1,i)).*(c2) +((r).^-1)'.*(c1.^2).*(w(:,1,i))) -B12(:,1,i).*((r).^-1)'.*(c1.^2).*(w(:,1,i)) -A55(:,1,i).*c1;
        f33=D11(:,1,i).*(c2+((r).^-1)'.*c1 -((r.^2).^-1)'.*c0) -A55(:,1,i).*c0;
        
        F_d([1:nr],[1:nr])=f11;
        F_d([1:nr],[nr+1:2*nr])=f12;
        F_d([1:nr],[2*nr+1:3*nr])=f13;
        F_d([nr+1:2*nr],[1:nr])=f21;
        F_d([nr+1:2*nr],[nr+1:2*nr])=f22;
        F_d([nr+1:2*nr],[2*nr+1:3*nr])=f23;
        F_d([2*nr+1:3*nr],[1:nr])=f31;
        F_d([2*nr+1:3*nr],[nr+1:2*nr])=f32;
        F_d([2*nr+1:3*nr],[2*nr+1:3*nr])=f33;
        
             % moshtaghat e sharte marzi for Circular plate with IM symply  (Circular  - IMS)
    
        F_d(1,:)=0; %u1=0
        F_d(1,1)=1;
        
        F_d(nr,:)=0;  %un=0;
        F_d(nr,nr)=1;
        
        F_d(nr+1,[1:nr])=(A11(1,1,i).*c1(1,:)).*(c1(1,:)*w(:,1,i)) +(A12(1,1,i).*c0(1,:)/b).*(c1(1,:)*w(:,1,i)); %%Qrz +Nrr*w,r  wrt w @r=0 ==0 
        F_d(nr+1,[nr+1:2*nr])=A55(1,1,i).*c1(1,:) +A11(1,1,i).*(c1(1,:)*u(:,1,i)).*c1(1,:) +A11(1,1,i).*(c1(1,:).*(c1(1,:)*w(:,1,i)).^2 +0.5*((c1(1,:)*w(:,1,i)).^2).*c1(1,:)) +A12(1,1,i).*(c0(1,:)*u(:,1,i)/b).*c1(1,:) +B11(1,1,i).*(c1(1,:)*phi(:,1,i)).*c1(1,:) +B12(1,1,i).*(c0(1,:)*phi(:,1,i)/b).*c1(1,:) -(c0(1,:)*N_TT(:,1,i))*c1(1,:); %Qrz +Nrr*w,r  wrt w @r=0 ==0 
        F_d(nr+1,[2*nr+1:3*nr])=A55(1,1,i).*c0(1,:) +B11(1,1,i).*c1(1,:).*(c1(1,:)*w(:,1,i)) +B12(1,1,i).*(c0(1,:)/b).*(c1(1,:)*w(:,1,i)); %Qrz +Nrr*w,r  wrt w @r=0 ==0 
        
        
        F_d(2*nr,:)=0;  %wn=0;
        F_d(2*nr,2*nr)=1;
        
        F_d(2*nr+1,:)=0;  %Phi1=0;
        F_d(2*nr+1,2*nr+1)=1;
        
        F_d(3*nr,[1:nr])=B11(nr,1,i).*c1(nr,:) +B12(nr,1,i).*c0(nr,:)/a; %M_rr wrt u @r=b ==0 
        F_d(3*nr,[nr+1:2*nr])=B11(nr,1,i).*(c1(nr,:)*w(:,1,i)).*c1(nr,:); %M_rr wrt w @r=b ==0 
        F_d(3*nr,[2*nr+1:3*nr])=D11(nr,1,i).*c1(nr,:) +D12(nr,1,i).*c0(nr,:)/a; %M_rr wrt phi @r=b ==0 
        
           deltaD=-inv(F_d)*FF_total(:,1);
           qq(:,1,i)=qq(:,1,i)+deltaD;
             
           sum1=0;
           sum2=0;
         
           for jj=1:length(deltaD)
               sum1=sum1+(deltaD(jj))^2;
               sum2=sum2+(qq(jj,1,i))^2;
               
           end
        
           itera=itera+1;
          
           
          u(:,1,i)=qq(1:nr,1,i);
          w(:,1,i)=qq(nr+1:2*nr,1,i);
          phi(:,1,i)=qq(2*nr+1:3*nr,1,i);
          
            if itera>2000
               break;
           end
        
   end
   error(i)=itera;
   
   QQ(:,1,i)=qq(1:3*nr,1,i);  % all deformation components
   WWW(:,1,i)=qq(nr+1:2*nr,1,i); 
   uuu(:,1,i)=qq(1:nr,1,i);
   u_mid(i,1)=uuu(ceil(nr/2),1,i);
   Wmid(i,1)=WWW(ceil(nr/2),1,i);
   Wfirst(i,1)=WWW(1,1,i);    
   u_first(i,1)=uuu(1,1,i);
       
    if (Wfirst(i,1)/h)>1
       
       break;
   end
   
  
   i
end 

figure(1)
%plot(tt,W1_linear/h)
hold on
plot(tt,W1_linear/h,'LineWidth', 2)
xlabel('time (sec)', 'Fontsize', 14, 'Color', 'black', 'FontName', 'Arial');
ylabel('W_{mid}/h', 'Fontsize', 14, 'Color', 'black', 'FontName', 'Arial');
title('W_{mid}/h vs time for different partial loading cases'); % Add a figure title
% Get the current axis and set the font size, color, and type
% ax = gca;
% set(ax,'Fontsize',14, 'color', 'black', 'FontName', 'Arial');
% dataset_name = sprintf('Load factor : %d%%', 100*r(length(index_heated_top)));
% legend({dataset_name}); 


figure(2)
hold on 
plot(r,reshape(W_linear(:,1,end),[nr,1,1])/h,'LineWidth', 2)
xlabel('r: Distance from the center of plate (m)', 'Fontsize', 14, 'Color', 'black', 'FontName', 'Arial');
ylabel('W/h', 'Fontsize', 14, 'Color', 'black', 'FontName', 'Arial');
title('W/h vs distance from the center for different partial loading cases'); % Add a figure title
% ax = gca;
% set(ax,'Fontsize',14, 'FontColor', 'black', 'FontName', 'Arial');
% dataset_name = sprintf('Load factor : %d%%', 100*r(length(index_heated_top)));
% legend({dataset_name}); 


[R,Z]=meshgrid(r,zz);
[N_rr,N_theta,M_rr,M_theta,sigma_rr,sigma_theta,sigma_shear,epsilon_rr,epsilon_theta,epsilon_shear,sigma_rr_mech,sigma_rr_thermal,sigma_theta_mech,sigma_theta_thermal] = post_processing_transient(nr,dt,a,b,zz,nz,h,k_shear,kisi,qq,right_edge,top_edge,bottom_edge,left_edge,partial_ratio_top,partial_ratio_bottom,A11,A22,A12,A55,B11,B22,B12,D11,D22,D12,N_TT,M_TT);
figure(4)
hold on 
plot(reshape(sigma_theta(:,1,end),[nz,1,1])/10^6,zz)
% hold on 
% plot(r,reshape(sigma_theta_mech(8,:,end),[1,nr,1])/10^6)
% hold on 
% plot(r,reshape(sigma_theta_thermal(8,:,end),[1,nr,1])/10^6)
% 
% figure(6)
% hold on 
% plot(r,reshape(sigma_rr(8,:,end),[1,nr,1])/10^6)
% hold on 
% plot(r,reshape(sigma_rr_mech(8,:,end),[1,nr,1])/10^6)
% hold on 
% plot(r,reshape(sigma_rr_thermal(8,:,end),[1,nr,1])/10^6)
% 
% figure(7)
% hold off 
% [M,c] = contour(R,Z,sigma_rr(:,:,end)/10^6,'ShowText','on');
% c.LineWidth = 3;
% 
% figure(8)
% hold off 
% [M,c] = contour(R,Z,sigma_theta(:,:,end)/10^6,'ShowText','on');
% c.LineWidth = 3;
