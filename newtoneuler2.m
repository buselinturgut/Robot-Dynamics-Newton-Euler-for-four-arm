clear all; clc;

syms the_i a_i alp_i d_i T(alp_i,a_i,d_i,the_i) l1 l2 l3 l4 l5 l6
syms t1 t2 t3 t4 t5 t6 t7 d1 d2 d3 d4 d5 d6 d7 h1 TT TP_t
syms r11 r12 r13 r21 r22 r23 r31 r32 r33 px py pz M
syms a b c d e f a1 a2 a3 a4
syms g
syms d_th1 d_th2 d_d3
syms dd_th1 dd_th2 dd_d3 m1 m2 m3
 
T(alp_i,a_i,d_i,the_i) = [cos(the_i)  -sin(the_i)    0   a_i;
    sin(the_i)*cos(alp_i)   cos(the_i)*cos(alp_i)   -sin(alp_i)   -sin(alp_i)*d_i;
    sin(the_i)*sin(alp_i)   cos(the_i)*sin(alp_i)   cos(alp_i)    cos(alp_i)*d_i;
    0                           0                         0           1];
    
    
D_H = [0 l1 0 t1;
       0 l2 0 t2];
 
dof = length(D_H(:,1));
TT = cell(dof,1);
TTS = TT;
T_inv = cell(dof,1);
TP = eye(4);
TP_st = cell(dof,1); %TP total yani aralarda daa ilgili eklme kadar ara ?ap?m bulma
T_tr = cell(dof,1);
%HP = eye(4);

for i = 1:dof
X = sprintf('transformation matrix %d - %d',i-1,i);
disp(X)
TT{i} = T(D_H(i,1),D_H(i,2),D_H(i,3),D_H(i,4));
TT{i}
TTS{i} = simplify(TT{i});
TP = TP * TT{i};
TP_st{i} = TP;
T_tr{i} = TTS{i}(1:3,1:3).';
T_inv{i} = [T_tr{i} -T_tr{i}*TT{i}(1:3,4);0 0 0 1];
end
disp('forward kinematics')
TP

w0 = [0 0 0]'; %0 an?ndaki a??sal h?z
dw0 = [0 0 0]'; %a??sal ivme
v0 = [0 0 0]'; %lineer h?z
dv0 =[0 -g 0].'; %yer?ekimi ivmesi
Z = [0 0 1]';

Dp = {d_th1 d_th2 d_d3}; %teta1dat, teta2dat..(birinci t?rev)
DDp = {dd_th1 dd_th2 dd_d3}; %datdat(ikinci t?rev)
Pc = [l1 0 0;...
     l2 0 0];
 
 M = [m1 m2 m3];
 Type = ['R', 'P', 'P'];

for i = 1:3
     syms(['Ixx',num2str(i)])
     syms(['Iyy',num2str(i)])
     syms(['Izz',num2str(i)])
end
 
for i = 1:3
     %3eklemlide ?u an i?in 2 de yaz?l?r
   % Ic{i} = [eval(['Ixx',num2str(i)])     0                        0;...
                %  0                       eval(['Iyy',num2str(i)]) 0;...
                %  0                       0                      eval(['Izz',num2str(i)])];
 
 %B?T?N elemanlar? 0 yapmak i?in 0 yap
 
 Ic{i} = [0 0 0;...
          0 0 0;...
          0 0 0];
 end
 
S =@(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

for i = 0:dof-1
    RT{i+1} = TT{i+1}(1:3,1:3).';
    P{i+1} = TT{i+1}(1:3,4);
  if strcmp(Type(i+1),'R')
    if i == 0
        w{i+1} = RT{i+1}*w0 + Dp{i+1}*Z;  %ilk turda bilgiler  burda kullan?l?yor 0 de?erinde yani(d?nel)
        dw{i+1} = RT{i+1}*dw0 + DDp{i+1}*Z+cross(RT{i+1}*w0,(Dp{i+1}*Z));
        dv{i+1} = RT{i+1}*(cross(dw0,P{i+1}) + cross(w0,cross(w0,P{i+1})) + dv0);
        dvc{i+1} = S(dw{i+1})*Pc(i+1,:).'+S(w{i+1})*S(w{i+1})*Pc(i+1,:).'+dv{i+1}; %k?tle merkezi do?rusal ivme
        F{i+1} = M(i+1)*dvc{i+1}; %K?tle merkezine ait kuvvet vekt?r?
        N{i+1} = Ic{i+1}*dw{i+1} + S(w{i+1})*(Ic{i+1}*w{i+1});
    else
        w{i+1} = RT{i+1}*w{i} + Dp{i+1}*Z;
        dw{i+1} = RT{i+1}*dw{i} + DDp{i+1}*Z+cross(RT{i+1}*w{i},(Dp{i+1}*Z));
        dv{i+1} = RT{i+1}*(cross(dw{i},P{i+1}) + cross(w{i},cross(w{i},P{i+1})) +dv{i});
        dvc{i+1} = S(dw{i+1})*Pc(i+1,:).'+S(w{i+1})*S(w{i+1})*Pc(i+1,:).'+dv{i+1};
        F{i+1} = M(i+1)*dvc{i+1};
        N{i+1} = Ic{i+1}*dw{i+1} + S(w{i+1})*(Ic{i+1}*w{i+1});
    end
  else
      if i == 0
        w{i+1} = RT{i+1}*w0;  %ilk turda bilgiler  burda kullan?l?yor 0 de?erinde yani(d?nel)
        dw{i+1} = RT{i+1}*dw0;
        dv{i+1} = RT{i+1}*(cross(dw{i},P{i+1}) + cross(w{i},cross(w{i},P{i+1})) +dv{i})+2*S(w{i+1})*(Dp{i+1}*Z)+DDp{i+1}*Z;
        dvc{i+1} = S(dw{i+1})*Pc(i+1,:).'+S(w{i+1})*S(w{i+1})*Pc(i+1,:).'+dv{i+1}; %k?tle merkezi do?rusal ivme
        F{i+1} = M(i+1)*dvc{i+1}; %K?tle merkezine ait kuvvet vekt?r?
        N{i+1} = Ic{i+1}*dw{i+1} + S(w{i+1})*(Ic{i+1}*w{i+1});
    else
        w{i+1} = RT{i+1}*w{i};
        dw{i+1} = RT{i+1}*dw{i};
        dv{i+1} = RT{i+1}*(cross(dw{i},P{i+1}) + cross(w{i},cross(w{i},P{i+1})) +dv{i})+2*S(w{i+1})*(Dp{i+1}*Z)+DDp{i+1}*Z;
        dvc{i+1} = S(dw{i+1})*Pc(i+1,:).'+S(w{i+1})*S(w{i+1})*Pc(i+1,:).'+dv{i+1};
        F{i+1} = M(i+1)*dvc{i+1};
        N{i+1} = Ic{i+1}*dw{i+1} + S(w{i+1})*(Ic{i+1}*w{i+1});
    end
      
      
      
      
  end
end
    
    syms fx fy fz nx ny nz  %i?e d?n?k denklemler
    MyRot = eye(3);
    %Myf = [fx;fy;fz];  %d??ardan etkileyen kuvvetler mylar
    Myf = [0;0;0];
    %Myn = [nx; ny; ny];
    Myn = [0; 0; 0];
    MyP = [0 0 0];
    
for i = dof:-1:1
        R{i} = TT{i}(1:3,1:3);
        if i == dof
           ff{i} = MyRot*Myf + F{i};
           nn{i} = N{i} + MyRot*Myn + S(Pc(i,:).')*F{i} + S(MyP)*MyRot*Myf;
        else 
            ff{i} = R{i+1}*ff{i+1} + F{i};
            nn{i} = N{i} + R{i+1}*nn{i+1} + S(Pc(i,:).')*F{i} + S(MyP)*R{i+1}*ff{i+1};
        end
end
    
    for i = 1:dof
        Tau{i} = nn{i}.'*Z;  %d?nel.  prizmatik ise fi*Z YAPILIR
    end