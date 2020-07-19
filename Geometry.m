%%LOAD IMAGE
im = imread('First.jpg');
baseImage = histeq(double(im)/255);

imshow(baseImage);

%[xi,yi] = getpts;
%disp(xi);
%disp(yi);
%% POINTS DEFINITION
hold on;

p1 = [546 0 1];
p2 = [1218 742 1];
p3 = [1589 481 1];
p4 = [2002 185 1];
p5 = [935 3029 1];
p6 = [1795 2580 1];
p7 = [2955 1968 1];
p8 = [11 2595 1];
p9 = [5 609 1];
p10 = [1125 1572 1];
p11 = [1649 1241 1];
p12 = [2387 761 1];
p13 = [2009 1637 1];
p14 = [2689 1241 1];
p16 = [2321 2975 1];
p17 = [3325 2507 1];
p18 = [1043 2219 1];
p19 = [3503 767 1];
p20 = [2615 95 1];
p21 = [3839 1307 1];
p22 = [3737 53 1];
p23 = [4015 1125 1];


x=[p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1) p8(1) p9(1) p10(1) p11(1) p12(1) p13(1) p14(1) p16(1) p17(1) p18(1) p19(1) p20(1) p21(1) p22(1) p23(1)];
y=[p1(2) p2(2) p3(2) p4(2) p5(2) p6(2) p7(2) p8(2) p9(2) p10(2) p11(2) p12(2) p13(2) p14(2) p16(2) p17(2) p18(2) p19(2) p20(2) p21(2) p22(2) p23(2)];
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); % plots points clicked by user with red circles

FNT_SZ = 16;
text(p1(1), p1(2), 'p1', 'FontSize', FNT_SZ, 'Color', 'w')
text(p2(1), p2(2), 'p2', 'FontSize', FNT_SZ, 'Color', 'w')
text(p3(1), p3(2), 'p3', 'FontSize', FNT_SZ, 'Color', 'w')
text(p4(1), p4(2), 'p4', 'FontSize', FNT_SZ, 'Color', 'w')
text(p5(1), p5(2), 'p5', 'FontSize', FNT_SZ, 'Color', 'w')
text(p6(1), p6(2), 'p6', 'FontSize', FNT_SZ, 'Color', 'w')
text(p7(1), p7(2), 'p7', 'FontSize', FNT_SZ, 'Color', 'w')
text(p8(1), p8(2), 'p8', 'FontSize', FNT_SZ, 'Color', 'w')
text(p9(1), p9(2), 'p9', 'FontSize', FNT_SZ, 'Color', 'w')
text(p10(1), p10(2), 'p10', 'FontSize', FNT_SZ, 'Color', 'w')
text(p11(1), p11(2), 'p11', 'FontSize', FNT_SZ, 'Color', 'w')
text(p12(1), p12(2), 'p12', 'FontSize', FNT_SZ, 'Color', 'w')
text(p13(1), p13(2), 'p13', 'FontSize', FNT_SZ, 'Color', 'w')
text(p14(1), p14(2), 'p14', 'FontSize', FNT_SZ, 'Color', 'w')
text(p16(1), p16(2), 'p16', 'FontSize', FNT_SZ, 'Color', 'w')
text(p17(1), p17(2), 'p17', 'FontSize', FNT_SZ, 'Color', 'w')
text(p18(1), p18(2), 'p18', 'FontSize', FNT_SZ, 'Color', 'w')
text(p19(1), p19(2), 'p19', 'FontSize', FNT_SZ, 'Color', 'w')
text(p20(1), p20(2), 'p20', 'FontSize', FNT_SZ, 'Color', 'w')
text(p21(1), p21(2), 'p21', 'FontSize', FNT_SZ, 'Color', 'w')
text(p22(1), p22(2), 'p22', 'FontSize', FNT_SZ, 'Color', 'w')
text(p23(1), p23(2), 'p23', 'FontSize', FNT_SZ, 'Color', 'w')

%% COMPUTE LINES FOR ALL 3 DIRECTIONS

%d1 lines (Right)
l1_1 = cross(p2, p4);
l1_2 = cross(p10, p12);
l1_3 = cross(p18, p19);
l1_4 = cross(p5, p7);
l1_5 = cross(p16, p17);

line([p2(1) p4(1)],[p2(2) p4(2)]);
line([p10(1) p12(1)],[p10(2) p12(2)]);
line([p18(1) p19(1)],[p18(2) p19(2)]);
line([p5(1) p7(1)],[p5(2) p7(2)]);
line([p16(1) p17(1)],[p16(2) p17(2)]);


%d2 lines (Forward)
l2_1 = cross(p1, p2);
l2_2 = cross(p9, p10);
l2_3 = cross(p8, p5);
l2_4 = cross(p22, p23);

line([p1(1) p2(1)],[p1(2) p2(2)]);
line([p9(1) p10(1)],[p9(2) p10(2)]);
line([p8(1) p5(1)],[p8(2) p5(2)]);
line([p22(1) p23(1)],[p22(2) p23(2)]);

%d3 lines (Up)
l3_1 = cross(p5, p2);
l3_2 = cross(p6, p3);
l3_3 = cross(p7, p4);
l3_4 = cross(p21, p20);

line([p5(1) p2(1)],[p5(2) p2(2)]);
line([p6(1) p3(1)],[p6(2) p3(2)]);
line([p7(1) p4(1)],[p7(2) p4(2)]);
line([p21(1) p20(1)],[p21(2) p20(2)]);

%% VANISHING POINTS
%d1
v1 = cross(l1_1,l1_2);
v1 = v1./v1(3);
v2 = cross(l1_2,l1_3);
v2 = v2./v2(3);
v3 = cross(l1_3,l1_4);
v3 = v3./v3(3);
v4 = cross(l1_4,l1_1);
v4 = v4./v4(3);
v5 = cross(l1_5,l1_1);
v5 = v5./v5(3);

vd1 = v1+v2+v3+v4+v5;
vd1 = vd1./5;

%d2
v1 = cross(l2_1,l2_2);
v1 = v1./v1(3);
v2 = cross(l2_2,l2_3);
v2 = v2./v2(3);
v3 = cross(l2_3,l2_4);
v3 = v3./v3(3);
v4 = cross(l2_4,l2_1);
v4 = v4./v4(3);

vd2 = v1+v2+v3+v4;
vd2 = vd2./4;

%d3
v1 = cross(l3_1,l3_2);
v1 = v1./v1(3);
v2 = cross(l3_2,l3_3);
v2 = v2./v2(3);
v3 = cross(l3_3,l3_4);
v3 = v3./v3(3);
v4 = cross(l3_4,l3_1);
v4 = v4./v4(3);

vd3 = v1+v2+v3+v4;
vd3 = vd3./4;

%Display vanishing points

x=[vd1(1) vd2(1) vd3(1)];
y=[vd1(2) vd2(2) vd3(2)];
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); % plots points clicked by user with red circles

FNT_SZ = 16;
text(vd1(1), vd1(2), 'vd1', 'FontSize', FNT_SZ, 'Color', 'w')
text(vd2(1), vd2(2), 'vd2', 'FontSize', FNT_SZ, 'Color', 'w')
text(vd3(1), vd3(2), 'vd3', 'FontSize', FNT_SZ, 'Color', 'w')

disp(vd1);
disp(vd2);
disp(vd3);


line([p7(1) vd1(1)],[p7(2) vd1(2)],'LineWidth',1.5);
line([p17(1) vd1(1)],[p17(2) vd1(2)],'LineWidth',1.5);
line([p19(1) vd1(1)],[p19(2) vd1(2)],'LineWidth',1.5);
line([p4(1) vd1(1)],[p4(2) vd1(2)],'LineWidth',1.5);
line([p12(1) vd1(1)],[p12(2) vd1(2)],'LineWidth',1.5);

line([p1(1) vd2(1)],[p1(2) vd2(2)],'LineWidth',1.5);
line([p9(1) vd2(1)],[p9(2) vd2(2)],'LineWidth',1.5);
line([p8(1) vd2(1)],[p8(2) vd2(2)],'LineWidth',1.5);
line([p6(1) vd2(1)],[p6(2) vd2(2)],'LineWidth',1.5);
line([p22(1) vd2(1)],[p22(2) vd2(2)],'LineWidth',1.5);

line([p5(1) vd3(1)],[p5(2) vd3(2)],'LineWidth',1.5);
line([p6(1) vd3(1)],[p6(2) vd3(2)],'LineWidth',1.5);
line([p7(1) vd3(1)],[p7(2) vd3(2)],'LineWidth',1.5);
line([p21(1) vd3(1)],[p21(2) vd3(2)],'LineWidth',1.5);

%% METRIC RECONSTRUCTION
% Definition of feature points on the plane than we want to metrically
% reconstruct

fp1 = [442 1959 1];
fp2 = [700 2125 1];
fp3 = [1052 2341 1];
fp4 = [1153 2283 1];
fp5 = [1383 2150 1];
fp6 = [3150 1127 1];
fp7 = [3477 938 1];
fp8 = [3609 861 1];
fp9 = [3475 430 1];
fp10 = [3360 70 1];

x=[fp1(1) fp2(1) fp3(1) fp4(1) fp5(1) fp6(1) fp7(1) fp8(1) fp9(1) fp10(1)];
y=[fp1(2) fp2(2) fp3(2) fp4(2) fp5(2) fp6(2) fp7(2) fp8(2) fp9(2) fp10(2)];

hold on;
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); 

FNT_SZ = 16;
text(fp1(1), fp1(2), 'fp1', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp2(1), fp2(2), 'fp2', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp3(1), fp3(2), 'fp3', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp4(1), fp4(2), 'fp4', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp5(1), fp5(2), 'fp5', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp6(1), fp6(2), 'fp6', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp7(1), fp7(2), 'fp7', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp8(1), fp8(2), 'fp8', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp9(1), fp9(2), 'fp9', 'FontSize', FNT_SZ, 'Color', 'w')
text(fp10(1), fp10(2), 'fp10', 'FontSize', FNT_SZ, 'Color', 'w')


line([fp1(1) fp3(1)],[fp1(2) fp3(2)]);
line([fp8(1) fp3(1)],[fp8(2) fp3(2)]);
line([fp8(1) fp10(1)],[fp8(2) fp10(2)]);
line([fp1(1) fp10(1)],[fp1(2) fp10(2)]);

r = norm([fp3(1), fp3(2)]-[fp8(1), fp8(2)]);

r1 = norm([fp3(1), fp3(2)]-[fp4(1), fp4(2)]) / r;
r2 = norm([fp3(1), fp3(2)]-[fp5(1), fp5(2)]) / r;
r3 = norm([fp3(1), fp3(2)]-[fp6(1), fp6(2)]) / r;
r4 = norm([fp3(1), fp3(2)]-[fp7(1), fp7(2)]) / r;

%Find corners of the two coplanar squares
sp11 = (1-r1) * fp2 + (r1) * fp9;
sp12 = (1-r2) * fp2 + (r2) * fp9;

sp13 = (1-r1) * fp1 + (r1) * fp10;
sp14 = (1-r2) * fp1 + (r2) * fp10;

sp21 = (1-r3) * fp2 + (r3) * fp9;
sp22 = (1-r4) * fp2 + (r4) * fp9;

sp23 = (1-r3) * fp1 + (r3) * fp10;
sp24 = (1-r4) * fp1 + (r4) * fp10;

x=[sp11(1) sp12(1) sp21(1) sp22(1) sp13(1) sp14(1) sp23(1) sp24(1)];
y=[sp11(2) sp12(2) sp21(2) sp22(2) sp13(2) sp14(2) sp23(2) sp24(2)];

hold on;
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); % plots points clicked by user with red circles

%% Affine rectification

linf = cross(vd1, vd2);
linf = linf./linf(3);

Hp1 = [1 0 0; 0 1 0; linf(1) linf(2) linf(3)];
affineImg = imwarp(baseImage,projective2d(Hp1'));
%%
imshow(affineImg);

%% Getting points of the two square after rectifying transformation
[x,y] = transformPointsForward(projective2d(Hp1'),[sp11(1) sp12(1) sp21(1) sp22(1) sp13(1) sp14(1) sp23(1) sp24(1)],[sp11(2) sp12(2) sp21(2) sp22(2) sp13(2) sp14(2) sp23(2) sp24(2)]);

sp11a = [x(1) y(1) 1];
sp12a = [x(2) y(2) 1];
sp13a = [x(5) y(5) 1];
sp14a = [x(6) y(6) 1];

sp21a = [x(3) y(3) 1];
sp22a = [x(4) y(4) 1];
sp23a = [x(7) y(7) 1];
sp24a = [x(8) y(8) 1];

hold on
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); % plots points clicked by user with red circles
%% Compute target positions for shape recontruction in the metrically reconstructed image
Np = 200;

sp13r = sp13a;
sp11r = [sp13r(1) sp13r(2)+Np 1];
sp12r = [sp13r(1)+Np sp13r(2)+Np 1];
sp14r = [sp13r(1)+Np sp13r(2) 1];

k = (sp13a - sp23a);
k = k([1 2]);
k = norm(k);

k1 = (sp13a - sp14a);
k1 = k1([1 2]);
k1 = norm(k1);

k = k/k1;

sp23r = [sp13r(1)+k*Np sp13r(2) 1];
sp21r = [sp23r(1) sp23r(2)+Np 1];
sp22r = [sp23r(1)+Np sp23r(2)+Np 1];
sp24r = [sp23r(1)+Np sp23r(2) 1];

x= [sp13r(1) sp11r(1) sp12r(1) sp14r(1) sp23r(1) sp21r(1) sp22r(1) sp24r(1)];
y= [sp13r(2) sp11r(2) sp12r(2) sp14r(2) sp23r(2) sp21r(2) sp22r(2) sp24r(2)];

hold on
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); 

%% Find shape reconstrution transformation with linear fit
syms a b c d e f;
Hr = [a b e; c d f; 0 0 1];

eq = [];
eq = [eq; sp11r' == Hr * sp11a'];
eq = [eq; sp12r' == Hr * sp12a'];
eq = [eq; sp13r' == Hr * sp13a'];
eq = [eq; sp14r' == Hr * sp14a'];

eq = [eq; sp21r' == Hr * sp21a'];
eq = [eq; sp22r' == Hr * sp22a'];
eq = [eq; sp23r' == Hr * sp23a'];
eq = [eq; sp24r' == Hr * sp24a'];

[A,y] = equationsToMatrix(eq,[a,b,c,d,e,f]);

X = [double(A)];
Y = [double(y)];

lm = fitlm(X,Y, 'y ~ x1 + x2 + x3 + x4 + x5 + x6 -1');
W = lm.Coefficients.Estimate;

Hr = [W(1,1) W(2,1) W(5,1); W(3,1) W(4,1) W(6,1); 0 0 1];
rectImg = imwarp(affineImg,projective2d(Hr.'));

%% Get final squares positions after shape reconstruction
imshow(rectImg);
[x,y] = transformPointsForward(projective2d(Hr.'),[sp11a(1) sp12a(1) sp21a(1) sp22a(1) sp13a(1) sp14a(1) sp23a(1) sp24a(1)],[sp11a(2) sp12a(2) sp21a(2) sp22a(2) sp13a(2) sp14a(2) sp23a(2) sp24a(2)]);
x = x+W(5,1);
y = y-W(6,1);
sp11r = [0 0 1];
sp12r = [0 0 1];
sp13r = [x(5) y(5) 1];
sp14r = [x(6) y(6) 1];

sp21r = [x(3) y(3) 1];
sp22r = [x(4) y(4) 1];
sp23r = [x(7) y(7) 1];
sp24r = [x(8) y(8) 1];

hold on
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); % plots points clicked by user with red circles


%% Get final feature point position in metrically encontructed plane and distance estimation
[x,y] = transformPointsForward(projective2d(Hp1.'),[fp1(1) fp2(1) fp3(1) fp4(1) fp5(1) fp6(1) fp7(1) fp8(1) fp9(1) fp10(1)],[fp1(2) fp2(2) fp3(2) fp4(2) fp5(2) fp6(2) fp7(2) fp8(2) fp9(2) fp10(2)]);
[x,y] = transformPointsForward(projective2d(Hr.'),x,y);
x = x+W(5,1);
y = y-W(6,1);
plot(x,y,'.r','MarkerSize',12, 'LineWidth', 3); 

fp1r = [x(1) y(1) 1];
fp2r = [x(2) y(2) 1];
fp3r = [x(3) y(3) 1];
fp4r = [x(4) y(4) 1];
fp5r = [x(5) y(5) 1];
fp6r = [x(6) y(6) 1];
fp7r = [x(7) y(7) 1];
fp8r = [x(8) y(8) 1];
fp9r = [x(9) y(9) 1];
fp10r = [x(10) y(10) 1];

disp(norm(fp1r-fp2r)/Np)
disp(norm(fp4r-fp5r)/Np)
disp(norm(fp6r-fp7r)/Np)
disp(norm(fp9r-fp10r)/Np)
disp(norm(fp3r-fp8r)/Np)

%% K MATRIX ESTIMATION
syms a b c d;
omega = [a 0 b; 0 1 c; b c d];

% Add constraints due to vanishing points
eqn = [];
eqn = [eqn; vd1*omega*vd2.' == 0];
eqn = [eqn; vd2*omega*vd3.' == 0];
eqn = [eqn; vd1*omega*vd3.' == 0];

% Add contraints due to homography
H = inv(Hr * Hp1);
h1 = H(:,1);
h2 = H(:,2);
eqn = [eqn; h1'*omega*h2 == 0];
eqn = [eqn; h1'*omega*h1 == h2'*omega*h2];

[A,y] = equationsToMatrix(eqn,[a,b,c,d]);

X = [double(A)];
Y = [double(y)];

lm = fitlm(X,Y, 'y ~ x1 + x2 + x3 + x4 - 1');

W = lm.Coefficients.Estimate;
IAC = double([W(1,1) 0 W(2,1); 0 1 W(3,1); W(2,1) W(3,1) W(4,1)]);

%Get parametrization values from IAC
alfa = sqrt(IAC(1,1));
u0 = -IAC(1,3)/(alfa^2);
v0 = -IAC(2,3);
fy = sqrt(IAC(3,3) - (alfa^2)*(u0^2) - (v0^2));
fx = fy /alfa;

% build K using the parametrization
K = [fx 0 u0; 0 fy v0; 0 0 1];

%% LOCALIZE CAMERA
x_ul = [0 0];
x_dl = [0 Np];
x_ur = [Np 0];
x_dr = [Np Np];

sp13ul = sp13;
sp11dl = sp11;
sp14ur = sp14;
sp12dr = sp12;

%Find transformation between world points and image points
H_omog = fitgeotrans([x_ul; x_dl; x_ur; x_dr], [sp13ul(1:2); sp11dl(1:2); sp14ur(1:2); sp12dr(1:2)], 'projective');
H_omog = H_omog.T.';

% Find rotation and translation matrix from H_omog
R = K \ H_omog;
lambda = 1 / norm(K \ H_omog(:,1));
R = R .* lambda;
R(:,3) = cross(R(:,1),R(:,2));

T = (K \ (lambda * h3));

% Find camera rotation and translation 
cameraRotation = R.';
cameraPosition = -R.'*T;
 
%% Vertical rectification
P = K * [R,T]; % projection matrix

H_vert_l_sr = inv([P(:,1), P(:,3), P(:,4)]);

alfa = 180;
R_y = [cos(deg2rad(alfa)) 0 sin(deg2rad(alfa)); ...
             0 1         0; ...
    -sin(deg2rad(alfa)) 0 cos(deg2rad(alfa))];

H_vert_l_sr = R_y * H_vert_l_sr ;

vertRecImg = imwarp(baseImage,projective2d(H_vert_l_sr.'));
imshow(vertRecImg);
