clc 
clear


%-----------------------input data-----------------------------------------
density = 922;
nu=0.658e-6;    %-----water 40 degree of centigrade------------------------
mu = density*nu;
% -------------------------------------------------------------------------
U=0.025;        %---m/s----   
Re=20;
H=Re*nu/U       % ---  Re= U H/nu---

L=11*H          % fully developed assumption------     
po = 1.0e+5;    %Outlet pressure
nx = 70;        % number of node in x direction--
dx = L/nx;

ny=50           % number of node in y direction--
%---------------------mesh generation--------------------------------------
ddy=1.0*H/ny

y=linspace(0,H,(ny+2))

dy= ddy*sin(pi/H*y)

r=0
for i=1:ny+2
    r=r+dy(i)
end

for i=1:ny+2
    dy(i)=dy(i)+(H-r)/(ny+2)
end

%--initial value for vertical velociy and related  matrice or coefficient--
v = zeros(nx+2,ny+1);   
vv = zeros(nx+2,ny+1);
aP_v = zeros(nx+2,ny+1);
Av=zeros(ny-1,ny-1);
Bv=zeros(ny-1,1);

%----initial value for horizontal velociy and relatedmatrice or coefficient
u = zeros(nx+1,ny+2)+U;  
uu = zeros(nx+1,ny+2)+U;
aP_u = zeros(nx+1,ny+2); 
A=zeros(ny,ny);
B=zeros(ny,1);

%-----------------------------initial value for pressure ------------------
p = zeros(nx+2,ny+2)+po;  
pc = zeros(nx+2,ny+2); 


%--------------------------under relaxation coefficient--------------------
u_under = 0.3;
v_under = 0.3;
p_under = 0.01;

%----------------------------convergence limitation------------------------
converg = 1.0;
simple_converg = 6e-08;
u_converg = 6e-08;
v_converg = 6e-08;
p_converg = 6e-08;


it=0;
%--------------------------SIMPLE algorithm--------------------------------
while converg>simple_converg
    it=it+1;
    converg = 0.0;
    
    uu = u;
    vv = v;
 %------------------- u-mommentum equation----------------------------------

uuu=zeros(ny,1);
uumom_converg=zeros(nx-1,ny);
umom_converg = 1.0;
while umom_converg > u_converg%____________________________________________
    umom_converg= 0.0;
    
    for I=3:nx+1%----------------------------------------------------------
        
        i=I-1;
        for J=2:ny+1    %--------------------------------------------------
              j=J-1;
              
            if(J==ny+1)%---------------------------------------------------
              
            Dw = mu/dx*dy(J);
            Ds = mu/dy(J)*dx;
            De = mu/dx*dy(J);
            Dn = mu/dy(J+1)*dx;
                
            fw = density*dy(J)*(u(i,J)+u(i-1,J))/2.0;
            fe = density*dy(J)*(u(i+1,J)+u(i,J))/2.0;
            fs = density*dx*(v(I,j)+v(I-1,j))/2.0;
            fn = density*dx*(v(I-1,j+1)+v(I,j+1))/2.0;
        
            pe_w = fw/Dw;
            pe_e = fe/De;
            pe_n = fn/Dn;
            pe_s = fs/Ds;
        
            aW=Dw*max(0,(1-.5*abs(pe_w)))+max(fw,0);
            aE=De*max(0,(1-.5*abs(pe_e)))+max(-fe,0);
            aS=Ds*max(0,(1-.5*abs(pe_s)))+max(fs,0);
            aN=Dn*max(0,(1-.5*abs(pe_n)))+max(-fn,0);
            aP = aW+aE+aS+aN+(fe-fw)*dy(J)+(fn-fs)*dx+mu*dx/dy(J);
       
            A(J-1,J-2)  = -aS;
            A(J-1,J-1)  = (aP/u_under+aN);
            B(J-1,1) = aE*uu(i+1,J) + aW*uu(i-1,J) + (p(I-1,J)-p(I,J))*dy(J) + (1-u_under)*aP/u_under*u(i,J);
            aP_u(i,J) = aP;
            elseif(J==2)   %-----------------------------------------------
            Dw = mu/dx*dy(J);
            Ds = mu/dy(J)*dx;
            De = mu/dx*dy(J);
            Dn = mu/dy(J+1)*dx;
                
            fw = density*dy(J)*(u(i,J)+u(i-1,J))/2.0;
            fe = density*dy(J)*(u(i+1,J)+u(i,J))/2.0;
            fs = density*dx*(v(I,j)+v(I-1,j))/2.0;
            fn = density*dx*(v(I-1,j+1)+v(I,j+1))/2.0;
        
            pe_w = fw/Dw;
            pe_e = fe/De;
            pe_n = fn/Dn;
            pe_s = fs/Ds;
        
            aW=Dw*max(0,(1-.5*abs(pe_w)))+max(fw,0);
            aE=De*max(0,(1-.5*abs(pe_e)))+max(-fe,0);
            aN=Dn*max(0,(1-.5*abs(pe_n)))+max(-fn,0);
            aS=Ds*max(0,(1-.5*abs(pe_s)))+max(fs,0);
            aP = aW+aN+aE+aS+(fe-fw)*dy(J)+(fn-fs)*dx+mu*dx/dy(J);
   
            A(J-1,J-1) = (aP/u_under+aS);
            A(J-1,J) = -aN;
            B(J-1,1) = aE*uu(i+1,J) + aW*uu(i-1,J) + (p(I-1,J)-p(I,J))*dy(J) + (1-u_under)*aP/u_under*u(i,J);
            aP_u(i,J) = aP;
            else%----------------------------------------------------------
            Dw = mu/dx*dy(J);
            Ds = mu/dy(J)*dx;
            De = mu/dx*dy(J);
            Dn = mu/dy(J+1)*dx;
                
            fw = density*dy(J)*(u(i,J)+u(i-1,J))/2.0;
            fe = density*dy(J)*(u(i+1,J)+u(i,J))/2.0;
            fs = density*dx*(v(I,j)+v(I-1,j))/2.0;
            fn = density*dx*(v(I-1,j+1)+v(I,j+1))/2.0;
        
            pe_w = fw/Dw;
            pe_e = fe/De;
            pe_n = fn/Dn;
            pe_s = fs/Ds;
        
            aW=Dw*max(0,(1-.5*abs(pe_w)))+max(fw,0);
            aE=De*max(0,(1-.5*abs(pe_e)))+max(-fe,0);
            aS=Ds*max(0,(1-.5*abs(pe_s)))+max(fs,0);
            aN=Dn*max(0,(1-.5*abs(pe_n)))+max(-fn,0);
            aP = aW+aN+aE+aS+(fe-fw)*dy(J)+(fn-fs)*dx;

            A(J-1,J-2) = -aS;
            A(J-1,J-1) = aP/u_under;
            A(J-1,J)  = -aN;
            B(J-1,1) = aE*uu(i+1,J) + aW*uu(i-1,J) + (p(I-1,J)-p(I,J))*dy(J) + (1-u_under)*aP/u_under*u(i,J);
            aP_u(i,J) = aP;
        end%---------------------------------------------------------------
        
        end%---------------------------------------------------------------
    
        uuu = Thomas(A,B);
        
    
        for jj=2:ny+1
            uumom_converg(i-1,jj-1) = abs(uu(i,jj)-uuu(jj-1));
            uu(i,jj) = uuu(jj-1);
        end
        
    end
    umom_converg=max(max(uumom_converg));
end%_______________________________________________________________________
%-------------------------- v-mommentum equation---------------------------


vvv=zeros(ny-1,1);

vvmom_converg=zeros(nx,ny-1);
vmom_converg = 1.0;

while vmom_converg>v_converg%______________________________________________
    
    vmom_converg = 0.0;
    for I=2:nx+1
        i=I-1;
        for J=3:ny+1
            
            j=J-1;
            if (j==2)
                
            Dw = mu/dx*dy(j);
            Ds = mu/dy(j)*dx;
            De = mu/dx*dy(j);
            Dn = mu/dy(j+1)*dx; 
                
            fw = density*dy(j)*(u(i,J)+u(i,J-1))/2.0;
            fe = density*dy(j)*(u(i+1,J)+u(i+1,J-1))/2.0;
           
            fs = density*dx*(v(I,j)+v(I,j-1))/2.0;
            fn = density*dx*(v(I,j+1)+v(I,j))/2.0;
            
            pe_w = fw/Dw;
            pe_e = fe/De;
            pe_n = fn/Dn;
            pe_s = fs/Ds;
        
            aW=Dw*max(0,(1-.5*abs(pe_w)))+max(fw,0);
            aE=De*max(0,(1-.5*abs(pe_e)))+max(-fe,0);
            aS=Ds*max(0,(1-.5*abs(pe_s)))+max(fs,0);
            aN=Dn*max(0,(1-.5*abs(pe_n)))+max(-fn,0);
            aP = aW+aN+aE+aS+(fe-fw)*dy(j)+(fn-fs)*dx;
                                        
            Av(j-1,j-1) = aP/v_under;
            Av(j-1,j) = -aN;
            Bv(j-1,1) = aE*vv(I+1,j) + aW*vv(I-1,j) + (p(I,J-1)-p(I,J))*dx + (1-v_under)*aP/v_under*v(I,j);
            aP_v(I,j) = aP;
            elseif(j==ny)
                
            Dw = mu/dx*dy(j);
            Ds = mu/dy(j)*dx;
            De = mu/dx*dy(j);
            Dn = mu/dy(j+1)*dx;
            
            fw = density*dy(j)*(u(i,J)+u(i,J-1))/2.0;
            fe = density*dy(j)*(u(i+1,J)+u(i+1,J-1))/2.0;
            fs = density*dx*(v(I,j)+v(I,j-1))/2.0;
            fn = density*dx*(v(I,j+1)+v(I,j))/2.0;
        
            pe_w = fw/Dw;
            pe_e = fe/De;
            pe_n = fn/Dn;
            pe_s = fs/Ds;
        
            aW=Dw*max(0,(1-.5*abs(pe_w)))+max(fw,0);
            aE=De*max(0,(1-.5*abs(pe_e)))+max(-fe,0);
            aS=Ds*max(0,(1-.5*abs(pe_s)))+max(fs,0);
            aN=Dn*max(0,(1-.5*abs(pe_n)))+max(-fn,0);
            aP = aW+aN+aE+aS+(fe-fw)*dy(j)+(fn-fs)*dx;
         
            Av(j-1,j-2)= -aS;
            Av(j-1,j-1) = aP/v_under;
            Bv(j-1,1) = aE*vv(I+1,j) + aW*vv(I-1,j) + (p(I,J-1)-p(I,J))*dx + (1-v_under)*aP/v_under*v(I,j);
            aP_v(I,j) = aP;
            else%----------------------------------------------------------
                
            Dw = mu/dx*dy(j);
            Ds = mu/dy(j)*dx;
            De = mu/dx*dy(j);
            Dn = mu/dy(j+1)*dx;
                
            fw = density*dy(j)*(u(i,J)+u(i,J-1))/2.0;
            fe = density*dy(j)*(u(i+1,J)+u(i+1,J-1))/2.0;
            fs = density*dx*(v(I,j)+v(I,j-1))/2.0;
            fn = density*dx*(v(I,j+1)+v(I,j))/2.0;
        
            pe_w = fw/Dw;
            pe_e = fe/De;
            pe_n = fn/Dn;
            pe_s = fs/Ds;
        
            aW=Dw*max(0,(1-.5*abs(pe_w)))+max(fw,0);
            aE=De*max(0,(1-.5*abs(pe_e)))+max(-fe,0);
            aS=Ds*max(0,(1-.5*abs(pe_s)))+max(fs,0);
            aN=Dn*max(0,(1-.5*abs(pe_n)))+max(-fn,0);
            aP = aW+aN+aE+aS+(fe-fw)*dy(j)+(fn-fs)*dx;
            
            Av(j-1,j-2) = -aS;
            Av(j-1,j-1) = aP/v_under;
            Av(j-1,j) = -aN;
            Bv(j-1,1) = aE*vv(I+1,j) + aW*vv(I-1,j) + (p(I,J-1)-p(I,J))*dx + (1-v_under)*aP/v_under*v(I,j);
            aP_v(I,j) = aP;
            end%-----------------------------------------------------------
        end
    
        vvv = Thomas(Av,Bv);
       
            for jj=2:ny
            vvmom_converg(I-1,jj-1) = abs(vv(I,jj)-vvv(jj-1));
            vv(I,jj) = vvv(jj-1);
        end
        
        end
    vmom_converg=max(max(vvmom_converg));
end%_______________________________________________________________________
%-------------------------------ppppppp------------------------------------

ppp_converg = zeros(nx,ny);
pp_converg = 1.0;

while pp_converg>p_converg
    
    pp_converg = 0.0;
    for i=2:nx+1
        for j=2:ny+1
           if (j==2)
               
            aS = 0.0;
            aN = density*dx*dx*v_under/aP_v(i,j);
         if i==nx+1
                aE = 0.0;
            else
                aE = density*dy(j)*dy(j)*u_under/aP_u(i,j);
            end
            if i==2
                aW = 0.0;
            else
                aW = density*dy(j)*dy(j)*u_under/aP_u(i-1,j);
            end
            aP = aE+aW+aS+aN;
        
            
        
            A(j-1,j-1) = aP/p_under;
            A(j-1,j) = -aN;
            B(j-1,1) = aE*pc(i+1,j) + aW*pc(i-1,j) + density*uu(i-1,j)*dy(j) - density*uu(i,j)*dy(j) + density*vv(i,j-1)*dx - density*vv(i,j)*dx;
        
           elseif(j==ny+1)
            
            aS = density*dx*dx*v_under/aP_v(i,j-1);
            aN = 0.0;
        if i==nx+1
                aE = 0.0;
            else
                aE = density*dy(j)*dy(j)*u_under/aP_u(i,j);
            end
            if i==2
                aW = 0.0;
            else
                aW = density*dy(j)*dy(j)*u_under/aP_u(i-1,j);
            end
            aP = aE+aW+aS+aN;
        
            
        
            A(j-1,j-2) = -aS;
            A(j-1,j-1) = aP/p_under;
            
            B(j-1,1) = aE*pc(i+1,j) + aW*pc(i-1,j) + density*uu(i-1,j)*dy(j) - density*uu(i,j)*dy(j) + density*vv(i,j-1)*dx - density*vv(i,j)*dx;
           else
        
            
            aS = density*dx*dx*v_under/aP_v(i,j-1);
            aN = density*dx*dx*v_under/aP_v(i,j);
        if i==nx+1
                aE = 0.0;
            else
                aE = density*dy(j)*dy(j)*u_under/aP_u(i,j);
            end
            if i==2
                aW = 0.0;
            else
                aW = density*dy(j)*dy(j)*u_under/aP_u(i-1,j);
            end
            aP = aE+aW+aS+aN;
        
            
        
            A(j-1,j-2) = -aS;
            A(j-1,j-1) = aP/p_under;
            A(j-1,j) = -aN;
            B(j-1,1) = aE*pc(i+1,j) + aW*pc(i-1,j) + density*uu(i-1,j)*dy(j) - density*uu(i,j)*dy(j) + density*vv(i,j-1)*dx - density*vv(i,j)*dx;
           end
        end
    
        ans = Thomas(A,B);
        
            for jj=2:ny+1
            ppp_converg(i-1,jj-1) = abs(pc(i,jj)-ans(jj-1));
            pc(i,jj) = ans(jj-1);
        end
    end
    pp_converg=max(max(ppp_converg));
end %______________________________________________________________________
% ------------------------- correction-------------------------------------

   %------- horizontal velocity correction------
for j=2:ny+1
    for i=2:nx
        uu(i,j) = uu(i,j) + dy(j)*u_under/aP_u(i,j)*(pc(i,j)-pc(i+1,j));
    end
end
    udiff=abs(u-uu);
    udiffcon=max(max(udiff));
    % ------ vertical velocity correction------
for j=2:ny
    for i=2:nx+1
        vv(i,j) = vv(i,j) + dx*v_under/aP_v(i,j)*(pc(i,j)-pc(i,j+1));
    end
end
    vdiff=abs(v-vv);
    vdiffcon=max(max(vdiff));
    
    %------pressure correction-------
    for j=2:ny+1
    for i=2:nx+1
        p(i,j) =p(i,j) + pc(i,j);
    end
end
   
    uu(1,:)=U;   %------------------- %Inlet bounday condition---------------
    vv(:,1) = 0;
    vv(1,:) = 0;%-------------------v Wall boundary condition---------------
    vv(:,ny+1)=0;    %
    uu(:,1)=-uu(:,2);   %-----------u Wall boundary condition---------------          
    uu(:,ny+2) =-uu(:,ny+1);  
    vv(nx+2,:)=vv(nx+1,:);  %------------%Outflow boundary condition-----------
    uu(nx+1,:)=uu(nx,:);  
    
    u = uu;
    v = vv;
    p(:,1) = p(:,2);       % topwall B.C
    p(1,:) = p(2,:);       %inlet pressure B.C
    p(:,ny+2) = p(:,ny+1); %bottomwall B.C
    
    converg = max(udiffcon,vdiffcon)
    converge(it)=converg;
    
end
% 
    for i=1:nx+1
for j=1:ny+1
 
         u(i,j)=(u(i,j)+u(i,j+1))/2  ;
%           
%          v(i,j)=(v(i,j)+v(i+1,j))/2  ;
%          p(i,j)=(p(i,j)+p(i,j+1)+p(i+1,j)+p(i+1,j+1))/4 ;
%         
    end
  u(i,ny+2)=(u(i,ny+1)+0)/2;
%          p(i,ny+2)=(p(i,ny+1)+Pout+p(i+1,ny+1)+Pout)/4;
% 
end

plot(u(nx-5,:),1:ny+2,'--.b')
hold on
plot(u(round(nx/6),:),1:ny+2,'-.r')
hold on
plot(u(round(2*nx/5),:),1:ny+2,'-*g')
hold on
plot(u(10,:),1:ny+2,'->y')
hold on
plot(u(2,:),1:ny+2,'-<r')
hold on

% 
% plot(u(nx,1:ny+2)/U,rrr(2:ny+3)/H,'--g')
% hold on
% plot(u(round(nx/2),1:ny+2)/U,rrr(2:ny+3)/H,'->b')
% hold on
% 
% plot(u(3,1:ny+2)/U,rrr(2:ny+3)/H,'-y')
% hold on
% plot(u(2,1:ny+2)/U,rrr(2:ny+3)/H,'-.r')
% hold on
% 
% 
% legend('outlet','middle section','initial section','inlet')
% xlabel('dimensionless velocity = u/U')
% ylabel('dimensionless height = y/H')
%   rr=zeros(ny+2)  
% rrr(1)=0
% for i=2:ny+3
%     rrr(i)=rrr(i-1)+dy(i-1)
% end
% 
% 
% contour(pppp)
% contour(uuuu)
% contour(vvvv)
% xlabel('horizontal direction node')
% ylabel('vertical direction node')
 plot(500:it,converge(499:it-1)) 
%  xlabel('iteration') 
%  ylabel('convergence')
%  plot((0:nx)/(nx),(u(:,(ny+2)/2))/U)
%  xlabel('dimensionless length = x/L')
%  ylabel('dimensionless velocity = u/U')
%  
%  plot((0:nx+1)/(nx+1),(v(:,(ny+2)/2)))
% xlabel('dimensionless length = x/L')
%  ylabel('vertical velocity')
