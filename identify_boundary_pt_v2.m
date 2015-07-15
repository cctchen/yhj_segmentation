function bp = identify_boundary_pt_v2(I2)
%% identify contour nodes on DXA image
scale = 0.2646e-3;

bnd_pts_file = 'sample_bp.dat';
key_pts_file = 'sample_kp.dat';


[junk,threshold] = edge(I2,'roberts');% return the threshold value
fudgefactor = 1.50; % enter a fudgefactor according to specific image that is dealing by this program
I4 = edge(I2,'roberts',threshold*fudgefactor);% detect the edge

se0 = strel('square',30);% enter a dilation value
I4 = imdilate(I4,se0); % dilate the image according to the dilation value entered above
bw = double(I4==1);
bw1 = flipud(bw);
bw2 = bwmorph(bw,'bridge');
bw3 = bwmorph(bw2,'close');
bw4 = bwmorph(bw3,'fill');
bw5 = bwmorph(bw4,'shrink');
bw6 = bwmorph(bw5,'thin');
bw7 = bwmorph(bw6,'spur');
bw8 = bwmorph(bw7,'erode',8);
bw9 = bwmorph(bw8,'thin',inf);% line 23 to 32 trim the dilated edge to thin contour

Segout =I2;
%Segout(bw9) = 255;

figure(10000);
imshow(I2),axis on; grid on; hold on
set(gcf,'outerposition',get(0,'screensize'));
%xlabel(strrep(case_id,'_',' '));

title({'Please select 3 points to obtian the circular femur head:1st point is on the left-down;';...
    '2nd point is on the right-up; 3rd point is on the top of head'},'FontSize',16);
[c1,c2] = ginput(3); % pick 3 points: the x coordinates of the three points are stored in array c1, the y coordinates of the three points are stored in array c2
plot(c1,c2,'ro');
%title({'Please select two (2) points to define the narrowest cross-section of the femoral neck';'1st point at lower neck; 2nd point at the upper neck'},'FontSize',16);
%[nx,ny] = ginput(2);
%plot(nx,ny,'r+-','LineWidth',2);
%title({'Please select two (2) points to define the inter-trochanteric line';'1st point at the middle of small trochanter; 2nd point at the top of great trochanter'},'FontSize',16);
%[tx,ty] = ginput(2);
%plot(tx,ty,'r+-','LineWidth',2);
%title({'Please select two (2) points to define the sub-trochanteric line';'1st point at the middle of small trochanter; 2nd point at the top of great trochanter'},'FontSize',16);
%[sx,sy] = ginput(2);
%plot(sx,sy,'r+-','LineWidth',2);
% title({'Please select one point to calculate the neck length'},'FontSize',16);
% [lx,ly] = ginput(1);
hold off;

[y1,x1]= find(bw9==1);% line 43 to 120 extract the coordinates of all the detected boundary points using subdivision method

B = 0.5*(c1(1)+max(x1));

N2 = find(x1<=B);
y11 = y1(N2);
x11 = x1(N2);

N3 = find(y11>=c2(1));
y12 = y1(N3);
x12 = x1(N3);
[yt1, m1, n1] = unique(-y12);
xt1 = x12(m1);

N4 = find(x1>=c1(2)&y1<=c2(1));
y2 = y1(N4);
x2 = x1(N4);
[xt2, m2, n2] = unique(x2);
yt2 = -y2(m2);

N5 = find(x1>B);
y3 = y1(N5);
x3 = x1(N5);

N6 = find (y3>c2(1));
y31 = y3(N6);
x31 = x3(N6);
[y32, m3, n3] = unique(-y31);
x32 = x31(m3);
[yt3, m4] = sort(y32,'descend');
xt3 = x32(m4);

if((c2(1)==c2(2))&(c2(2)==c2(3)))
    disp('The selected points cannot form a circle');
elseif((c2(1)~=c2(2))&(c2(2)~=c2(3)))
    k1=(c1(2)-c1(1))/(c2(2)-c2(1));
    k2=(c1(3)-c1(2))/(c2(3)-c2(2));
end
if(k1==k2)
    disp('The selected points cannot form a circle');
end
a=2*(c1(2)-c1(1));
b=2*(c2(2)-c2(1));
c=c1(2)*c1(2)+c2(2)*c2(2)-c1(1)*c1(1)-c2(1)*c2(1);
d=2*(c1(3)-c1(2));
e=2*(c2(3)-c2(2));
f=c1(3)*c1(3)+c2(3)*c2(3)-c1(2)*c1(2)-c2(2)*c2(2);
%disp('Center:');
Rx=(b*f-e*c)/(b*d-e*a);
Ry=(d*c-a*f)/(b*d-e*a);
%disp('Radius');
R=sqrt((Rx-c1(1))*(Rx-c1(1))+(Ry-c2(1))*(Ry-c2(1)));

pc=[Rx Ry];
p1=[c1(1) c2(1)];
p2=[c1(2) c2(2)];
p3=[Rx+R Ry];
px1=p1-pc;
px2=p2-pc;
px3=p3-pc;
Ar=dot(px2,px1)/(norm(px1)*norm(px2));
Ar=acos(Ar);
Ars=dot(px1,px3)/(norm(px1)*norm(px3));
Ars=acos(Ars);
if Ar < pi
    Ar= 2*pi-Ar;
%else
end
seta=Ars:0.05:Ar+Ars;
Rxx = Rx+R*cos(seta);
Rxx = Rxx';
Ryy = Ry+R*sin(seta);
Ryy = Ryy';
[bo,to]=size(I2);
yt=[-bo;yt1;-Ryy;yt2;yt3;-bo];
xt=[xt1(1);xt1;Rxx;xt2;xt3;xt3(max(size(xt3)))];

%% ========================================================================
%% (1) find the narrowest cross-section at the neck
%% (2) find an intertrochanteric cross-section that is parallel to the
%% narrowest neck cross-section
%% (3) find a subtrochnateric cross-section

cx1 = c1(1); cy1 = -c2(1);
cx2 = c1(2); cy2 = -c2(2);
dd1 = sqrt((cx1 - xt).^2 + (cy1 - yt).^2);
dd2 = sqrt((cx2 - xt).^2 + (cy2 - yt).^2);
nnd = length(xt);

[tmp1,i1] = min(dd1);
[tmp2,i2] = min(dd2);

curve1 = [xt(1:i1) yt(1:i1)];       nn1 = length(curve1(:,1));
curve2 = [xt(i2:nnd) yt(i2:nnd)];   nn2 = length(curve2(:,1));

nk1 = [xt(i1) yt(i1)];
nk2 = [xt(i2) yt(i2)];
while (1) 
    dp1p2 = norm(nk1-nk2);
    dp1c2 = sqrt((nk1(1).*ones(nn2,1) - curve2(:,1)).^2 + (nk1(2).*ones(nn2,1) - curve2(:,2)).^2);
    [tmp, i2] = min(dp1c2);
    tmp_p2 = curve2(i2,:);
    tmp_dd = norm(nk1-tmp_p2);
    if tmp_dd < dp1p2
        dp1p2 = tmp_dd;
        nk2 = tmp_p2;
    else
        break;
    end;
    
    dp2c1 = sqrt((nk2(1).*ones(nn1,1) - curve1(:,1)).^2 + (nk2(2).*ones(nn1,1) - curve1(:,2)).^2);
    [tmp, i1] = min(dp2c1);
    tmp_p1 = curve1(i1,:);
    tmp_dd = norm(nk2-tmp_p1);
    if tmp_dd < dp1p2
        dp1p2 = tmp_dd;
        nk1 = tmp_p1;
    else
        break;
    end;
end;

wd = min([max(xt)-min(xt) max(yt)-min(yt)]);
dx = wd/4*abs(nk2(1)-nk1(1))/norm(nk2-nk1);
dy = wd/4*abs(nk2(2)-nk1(2))/norm(nk2-nk1);

tx1 = nk1(1) + dx; ty1 = nk1(2) - dy;
tx2 = nk2(1) + dx; ty2 = nk2(2) - dy;

tx1 = tx1 - dx; ty1 = ty1 - dy;
tx2 = tx2 + 2*dx; ty2 = ty2 + 2*dy;

nx = [nk1(1); nk2(1)]; ny = [nk1(2); nk2(2)];
tx = [tx1; tx2]; ty = [ty1; ty2];
sx = [tx1; tx2]; sy = [ty1; ty1];
%% ------------------------------------------------------------------------

xx1 = xt;
yy1 = yt;
xx1(4:2:length(xx1))=[];
yy1(4:2:length(yy1))=[];% line 122 to line 129 reduce the number of detected keypoints that will be inputted into ANSYS
xx1(2:2:length(xx1))=[];
yy1(2:2:length(yy1))=[];

yy1(length(yy1))=yy1(1);
[Xs Ys]=smooth_contours(xx1,yy1,5);
%[OXs OYs]=smooth_contours(xt,yt,1);
bp = [Xs Ys].*scale;
cp = [c1 -c2; nx ny; tx ty; sx sy].*scale;

%% make the femur bone in the 'after' and 'before' cases have the same
%% length
coord_trans = 0;
if coord_trans == 1
    if length(strfind(case_id,'after'))>0
        case_before_1 = [output_dir '\' strrep(case_id,'after','before') '_1.dat'];
        case_before_2 = [output_dir '\' strrep(case_id,'after','before') '_2.dat'];
        bp_bf = load(case_before_1);
        cp_bf = load(case_before_2);
        bp_af = bp;
        cp_af = cp;
        
        %% coordinate transformation parameters
        theta_bf = 0.; dx_bf = 0.; dy_bf = 0.;
        theta_af = 0.; dx_af = 0.; dy_af = 0.;
        
        %% do the rotation first here
        rotation = 0;  %% there is a bug in the rotation codes
        if rotation == 1
            v_xx = [1. 0.];
            v_bf = [cp_bf(5,1)-cp_bf(4,1) cp_bf(5,2)-cp_bf(4,2)];
            v_af = [cp_af(5,1)-cp_af(4,1) cp_af(5,2)-cp_af(4,2)];
            bf_xx = acos(dot(v_bf,v_xx)/(norm(v_bf)*norm(v_xx)));
            af_xx = acos(dot(v_af,v_xx)/(norm(v_af)*norm(v_xx)));
            theta = acos(dot(v_af,v_bf)/(norm(v_af)*norm(v_bf)));
            %theta = abs(bf_xx - af_xx);
            if (af_xx < bf_xx)
                bp_af = [bp_af(:,1).*cos(theta)-bp_af(:,2).*sin(theta) ...
                    bp_af(:,1).*sin(theta)+bp_af(:,2).*cos(theta)];
                cp_af = [cp_af(:,1).*cos(theta)-cp_af(:,2).*sin(theta) ...
                    cp_af(:,1).*sin(theta)+cp_af(:,2).*cos(theta)];
                theta_bf = 0.;
                theta_af = -theta;
            elseif (bf_xx < af_xx)
                bp_bf = [bp_bf(:,1).*cos(theta)-bp_bf(:,2).*sin(theta) ...
                    bp_bf(:,1).*sin(theta)+bp_bf(:,2).*cos(theta)];
                cp_bf = [cp_bf(:,1).*cos(theta)-cp_bf(:,2).*sin(theta) ...
                    cp_bf(:,1).*sin(theta)+cp_bf(:,2).*cos(theta)];
                theta_bf = -theta;
                theta_af = 0.;
            end;
        end;

        %% do the translation
        translation = 1;
        if translation ==1 
            dx = abs(min(bp_bf(:,1)) - min(bp_af(:,1)));
            dy = abs(max(bp_bf(:,2)) - max(bp_af(:,2)));
            if (min(bp_bf(:,1)) < min(bp_af(:,1)))
                bp_bf(:,1) = bp_bf(:,1) + dx;
                cp_bf(:,1) = cp_bf(:,1) + dx;
                dx_bf = -dx;
                dx_af = 0.;
            elseif (min(bp_af(:,1)) < min(bp_bf(:,1)))
                bp_af(:,1) = bp_af(:,1) + dx;
                cp_af(:,1) = cp_af(:,1) + dx;
                dx_bf = 0.;
                dx_af = -dx;
            end
            if (max(bp_bf(:,2)) < max(bp_af(:,2)))
                bp_bf(:,2) = bp_bf(:,2) + dy;
                cp_bf(:,2) = cp_bf(:,2) + dy;
                dy_bf = -dy;
                dy_af = 0.;
            elseif (max(bp_af(:,2)) < max(bp_bf(:,2)))
                bp_af(:,2) = bp_af(:,2) + dy;
                cp_af(:,2) = cp_af(:,2) + dy;
                dy_bf = 0.;
                dy_af = -dy;
            end;
            cp_bf(10:11,:) = [theta_bf 0.; dx_bf dy_bf];
        end;
        
        %% do the tailoring so the 'before' and 'after' have the same
        %% length
        tailoring = 1;
        if tailoring == 1
            h_bf = abs(max(bp_bf(:,2)) - min(bp_bf(:,2)));
            h_af = abs(max(bp_af(:,2)) - min(bp_af(:,2)));
            hh = min([h_bf h_af]);

            ind_bf = find(bp_bf(:,2)<(max(bp_bf(:,2))-hh));
            bp_bf(ind_bf,:) = [];
            ind_af = find(bp_af(:,2)<(max(bp_af(:,2))-hh));
            bp_af(ind_af,:) = [];
            
            bottom = max([min(bp_bf(:,2)) min(bp_af(:,2))]);
            %bottom = max([bp_bf(1,2) bp_bf(length(bp_bf(:,1)),2)...
            %              bp_af(1,2) bp_af(length(bp_af(:,1)),2)]);
            bp_bf(1,2) = bottom; bp_bf(length(bp_bf(:,1)),2) = bottom;
            bp_af(1,2) = bottom; bp_af(length(bp_af(:,1)),2) = bottom;

            bp = bp_af;
            cp(10:11,:) = [theta_af 0.; dx_af dy_af];

            if exist(case_before_1) ~= 0
                delete(case_before_1);
            end;
            dlmwrite(case_before_1,bp_bf,'delimiter', '\t','precision',10);

            if exist(case_before_2) ~= 0
                delete(case_before_2);
            end;
            
            dlmwrite(case_before_2,cp_bf,'delimiter', '\t','precision',10);
        end;
        if test_it==1
            figure(1000+case_i);
            plot(bp_bf(:,1),bp_bf(:,2),'b-.',bp_af(:,1),bp_af(:,2),'r-.'); axis equal; grid on;
        end;
    end;
end;

%dlmwrite(bnd_pts_file,bp,'delimiter', '\t','precision',10);
%dlmwrite(key_pts_file,cp,'delimiter', '\t','precision',10);
    