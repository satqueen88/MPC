
clear
MAXTIME = 3000;
x1 = [0;0;0;0;0]; %��ԕϐ�[��,��,x,y,��]�̏����l
xd = [0;4;0];%�ڕW�lgannma,y,thete
N = 10; %N���܂ŏ�ԗ\��
Dt=0.01;%��ԗ\���̊Ԋu
x_next=x1;
x_estimat=zeros(5,2);


for count = 1:MAXTIME
    %%optimization
    U = ones(N,1);
    %U = [5;5;5;5;5;5;5;5;5;5]
    %hatx=zeros(5,N+1);
    %hatx(:,1)=x1;
    
    if x_next(3,:)< 0
    xd(2)=0;
    elseif x_next(3,:)< 200
    xd(2)=x_next(3,:)*(2/100);
    else
    xd(2)=4;
    end
    
    
    
    fun = @(U)costfunction(x_next,xd,U,N,Dt,count);%���݂̈ʒu�ŕ]���֐��̌v�Z
    [upot,fval]=lsqnonlin(fun,U);
    [upot]=lsqnonlin(fun,U);%�]���֐����ŏ���������Ǌp(delta)��T��
    
    hatx = dynamics(x_next, upot, N,Dt,count);%�ŏ����������Ǌp�ő��ǂ��s�����ꍇ�̏�ԗ\�� 
    u(count,1)=upot(1);
    
    %plot(hatx(3,:),hatx(4,:),'-o');
	%title(['estimated trajectory']);
	%xlabel('x [m]');
	%ylabel('y [m]');
	%hold on;
    %grid on;
    %if count < 1000
    %x_estimat=missdynamics(hatx(:,1),upot(1),1,Dt,count);%���̈���(��,��,x,y,��)�̏�ԁA���̈���̎��Ǌp(delta)
    %else
    x_estimat=dynamics(hatx(:,1),upot(1),1,Dt,count)
    %end
    x_next=x_estimat(:,2);%���ۂɑ��ǂ��s���Ă݂�
    xy(count,1)=x_next(3,:);
    xy(count,2)=x_next(4,:);
    
    %if rem(count,100)==0
    plot(x_next(3,:),x_next(4,:),'-o');
    title(['observed trajectory']);
    xlabel('x [m]');
    ylabel('y [m]');
    hold on;
    %grid on;
    %end
end


function hatx = dynamics(x1,delta,N,Dt,count)
    V  = 100*1000/3600;%���������x
    %Dt = 0.001; %1ms
    m = 1500;%�ԗ��d��
    Kf = 55000;%�O�փR�[�i�����O�p���[
    Kr = 60000;%��փR�[�i�����O�p���[
    lf = 1.1;%�ԗ��d�S����O�ւ܂�
    lr = 1.6;%�ԗ��d�S�����ւ܂�
	I = 2500;%���[�C���O�������[�����g
       
    
    %m = 1500;
    %Kf = 55000;
    %Kr = 60000;
    %lf = 1.1;
    %lr = 1.6;
	%I = 2500;
    
    
    
    
	a11 = -2*(Kf+Kr)/(m*V);
	a12 = -(m*V+2/V*(lf*Kf-lr*Kr))/(m*V);
	a21 = -2*(lf*Kf-lr*Kr)/I;
	a22 = -2*(lf*lf*Kf+lr*lr*Kr)/(I*V);
	
	A = [a11, a12; a21, a22];
	B = [2*Kf/(m*V), 2*lf*Kf/I]';
    
    %hat = zeros(5,N+1);
    hatx = zeros(5,N+1);
    %hatx(:,1) = hatx1;
    hatx(:,1)=x1;
    for i = 1:N
        fu=zeros(5,1);
        fu(1) = A(1)*hatx(1,i) + B(1)*delta(i);%d/dt(beta)
        fu(2) = A(2).*hatx(2,i) + B(2)*delta(i);%d/dt(gammma)
        fu(3) = V*cos(hatx(5,i));%d/dt(x)
        fu(4) = V*sin(hatx(5,i));%d/dt(y)
        fu(5) = hatx(2,i);	%d/dt(theta)
        hatx(:,i+1) = hatx(:,i)  + fu*Dt;%10���܂ŗ\��
    end
end

function mhatx = missdynamics(x1,delta,N,Dt,count)
    V  = 30*1000/3600;%���������x
    %Dt = 0.001; %1ms
    m = 1500;%�ԗ��d��
    Kf = 55000;%�O�փR�[�i�����O�p���[
    Kr = 60000;%��փR�[�i�����O�p���[
    lf = 1.1;%�ԗ��d�S����O�ւ܂�
    lr = 1.6;%�ԗ��d�S�����ւ܂�
	I = 2500;%���[�C���O�������[�����g
       
    
    %m = 1500;
    %Kf = 55000;
    %Kr = 60000;
    %lf = 1.1;
    %lr = 1.6;
	%I = 2500;
    
    
    
    
	a11 = -2*(Kf+Kr)/(m*V);
	a12 = -(m*V+2/V*(lf*Kf-lr*Kr))/(m*V);
	a21 = -2*(lf*Kf-lr*Kr)/I;
	a22 = -2*(lf*lf*Kf+lr*lr*Kr)/(I*V);
	
	A = [a11, a12; a21, a22];
	B = [2*Kf/(m*V), 2*lf*Kf/I]';
    
    %hat = zeros(5,N+1);
    mhatx = zeros(5,N+1);
    %hatx(:,1) = hatx1;
    mhatx(:,1)=x1;
    for i = 1:N
        fu=zeros(5,1);
        fu(1) = A(1)*mhatx(1,i) + B(1)*delta(i);%d/dt(beta)
        fu(2) = A(2).*mhatx(2,i) + B(2)*delta(i);%d/dt(gammma)
        fu(3) = V*cos(mhatx(5,i));%d/dt(x)
        fu(4) = V*sin(mhatx(5,i));%d/dt(y)
        fu(5) = mhatx(2,i);	%d/dt(theta)
        mhatx(:,i+1) = mhatx(:,i)  + fu*Dt;%10���܂ŗ\��
    end
end


function	fvec = costfunction(x1, xd, umat, N, Dt, count)
%hat=zeros(5,N+1);
%if count < 100
%hat=missdynamics(x1, umat, N, Dt ,count); 
%else
hat=dynamics(x1, umat, N, Dt ,count); 
%end
%gamma_mat=zeros(1,N);
%y_mat=zeros(1,N);
%theta_mat=zeros(1:N);

%emat=zeros(3,N);
gamma_mat = xd(1)*ones(1,N) - hat(2,2:N+1);
y_mat = xd(2)*ones(1,N) - hat(4,2:N+1);
theta_mat = xd(3)*ones(1,N) - hat(5,2:N+1);
%emat(1,1:N) = xd(1)*ones(1,N) - hat(2,2:N+1);
%emat(2,1:N) = xd(2)*ones(1,N) - hat(4,2:N+1);
%emat(3,1:N) = xd(3)*ones(1,N) - hat(5,2:N+1);
%gamma_mat=emat(1);
%y_mat=emat(2);
%theta_mat=emat(3);

gamma_e=reshape(gamma_mat,N,1);%n,1�Ɍ`��ύX
y_e=reshape(y_mat,N,1);%n,1�Ɍ`��ύX
theta_e=reshape(theta_mat,N,1);%n,1�Ɍ`��ύX
%q=[10;10;10];
%qq=repmat(q,N,1);
%Q = diag(qq);
%e=reshape(emat,3*N,1);

a=10;%gammma_geinn
aa=repmat(a,N,1);%a��N,1�ɃR�s�[
A=diag(aa);%�Ίp�s��쐬
AA=chol(A);%�R���X�L�[���q�v�Z

b=10;%y_geiinn
bb=repmat(b,N,1);%a��N,1�ɃR�s�[
B=diag(bb);%�Ίp�s��쐬
BB=chol(B);%�R���X�L�[���q�v�Z

c=300;%theta_geinn
cc=repmat(c,N,1);%a��N,1�ɃR�s�[
C=diag(cc);%�Ίp�s��쐬
CC=chol(C);%�R���X�L�[���q�v�Z

d=200;%delta_geinn
dd=repmat(d,N-1,1);%a��N,1�ɃR�s�[
D = diag(dd);%�Ίp�s��쐬
DD=chol(D);%�R���X�L�[���q�v�Z

for i = 1:N-1
    du(i)=umat(i)-umat(i+1);%���͂̍�
end

%QQ=chol(Q);


%fvec = [QQ*e; R*du'];
fvec = [AA*gamma_e;BB*y_e;CC*theta_e; DD*du'];%�]���֐�J
end