function flag = is_overlap(t1,t2,r1,r2,H)
% IS_OVERLAP computes whether two conic overlap, if overlap, return 1,
% otherwise, return 0. The first conic is an ellipse which is a projection
% of an circle with centre T1 and  radius R1 by a homography H. The other
% conic is an circle with centre T2 and  radius R2.
% Then the two conic can be described as follow:
%                       [x]                                           [     1       0                   -t1(1)           ]                       [     1       0                   -t2(1)           ]
% [x, y, 1] * Ci * [y] = 0,     where C1=invHT* [     0       1                   -t1(2)           ] * invH,    C2=[     0       1                   -t2(2)           ].
%                       [1]                                           [-t1(1)  -t1(2)  t1(1)^2+t1(2)^2-r1^2]                       [-t2(1)  -t2(2)  t2(1)^2+t2(2)^2-r2^2]
%   We decide an pair of conics are overlap by first check whether one
%   conic centre is inside the other conic, if so, there overlap, return;
%   otherwise, we compute if the two conic intersect in succession.
%   If the equation set composed by the two conics have real root,
%   then they overlap.
%   |A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
%   |(x-t2(1))^2 + (y-t2(2))^2 = r2^2
%   where A=C1(1,1), B=C1(1,2)+C1(2,1), C=C1(2,2), D=C1(1,3)+C1(3,1),
%       E=C1(2,3)+C1(3,2), F=C1(3,3).
%   The equation set can be written as
%   |A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
%   |x^2 + y^2 = r2^2
%   where newA=A, newB=B, newC=C, newD=2*t2(1)*A+t2(2)*B+D,
%       newE=t2(1)*B+2*t2(2)*C+E, newF=t2(1)^2*A+t2(1)*t2(2)*B+t2(2)^2*C+t2(1)*D+t2(2)*E+F.
%   Then we get the quartic equation and check whether it have a real root,
%   c4*y^4+c3*y^3+c2*y^2+c1*y+c0=0
%   where c4=(C-A)^2+B^2, c3=2*(C*E+B*D-A*E), c2=E^2+D^2-r2^2*B^2+2*(C-A)*(A*r2^2+F),
%       c1=2*E*(A*r2^2+F)-2*B*D*r2^2, c0=(A*r2^2+F)^2-D^2*r2^2.

invH = inv(H);
flag = 0;
% Check if the centre of conic1 is inside the conic2
t = H*[t1(1); t1(2); 1];
t(1)=t(1)/t(3); t(2)=t(2)/t(3);
if (t(1)-t2(1))^2+(t(2)-t2(2))^2<r2^2
    flag = 1;
    return;
end
C1=[1 0 -t1(1); 0 1 -t1(2); -t1(1) -t1(2), t1(1)^2+t1(2)^2-r1^2];
C1_=invH'*C1*invH;
% Check if the centre of conic1 is inside the conic2
if([t2(1) t2(2) 1]*C1_*[t2(1); t2(2); 1])<0
    flag = 1;
    return;
end
A=C1_(1,1);
B=C1_(1,2)+C1_(2,1);
C=C1_(2,2);
D=C1_(1,3)+C1_(3,1);
E=C1_(2,3)+C1_(3,2);
F=C1_(3,3);
D_=2*t2(1)*A+t2(2)*B+D;
E_=t2(1)*B+2*t2(2)*C+E;
F_=t2(1)^2*A+t2(1)*t2(2)*B+t2(2)^2*C+t2(1)*D+t2(2)*E+F;
D=D_;
E=E_;
F=F_;
c4=(C-A)^2+B^2;
c3=2*(C*E-A*E+B*D);
c2=E^2+D^2-r2^2*B^2+2*(C-A)*(A*r2^2+F);
c1=2*E*(A*r2^2+F)-2*B*D*r2^2;
c0=(A*r2^2+F)^2-D^2*r2^2;
r = roots([c4 c3 c2 c1 c0]);
for i =1:length(r)
    if(isreal(r(i)))
        flag = 1;
        break;
    end
end
