clear; clc; close all;
global case_i output_dir scale test_it;
FUC=1.073;
INTERCEPT=-1024.3;
SLOP=1.2882;
loca=0;
Pthreshold=50;
mu=0.3;
pixel_A = (0.48828*0.001)^2;
pixel_scale=0.48828*0.001;

par = 2; % 1 -- calculate CT cross-section properties; 2 -- calculate CTXA cross-section properties
%CSA=9.892; %cm2;
path = 'E:\QCT_CS_data';
[TMP,CSV_NAME]=xlsread([path '\hipdb_0709.csv'],'D:D');
CSV_CSA=xlsread([path '\hipdb_0709.csv'],'FP:FP');
CSMI_inplane=xlsread([path '\hipdb_0709.csv'],'FX:FX');
CSMI_outplane = xlsread([path '\hipdb_0709.csv'],'FY:FY');
SM_inplane=xlsread([path '\hipdb_0709.csv'],'GA:GA');
if par == 1
    [status,list]=system('dir E:\QCT_CS_data\*FN.dcm /S/B');
elseif par == 2
    [status,list]=system('dir E:\QCT_CS_data\*CTXA.dcm /S/B');
end;

filelist = strsplit(list);
filelist(strcmp(filelist,''))=[];

[file_temp,filenumber]=size(filelist);

if par == 1
    EA(filenumber)=0;
    GA(filenumber)=0;
    EIx(filenumber)=0;
    EIy(filenumber)=0;
    
    for FN=1:filenumber
        pixel = dicomread(filelist{FN});
        [mm,nn] = find(pixel>1);
        pixel([mm nn]) = pixel([mm nn]) + 1024; 
        pixel(:,1) = 0;
        
        rou=(double(FUC).*(double(pixel)))./double(SLOP);
        E=double(29.8).*(double((1e-3)).*double(rou)).^1.56.*1.0e9;
        
        EA(FN)=sum(sum(double(E).*pixel_A));
        GA(FN)=EA(FN)./(2*(1+mu));
        
        x_c = (max(nn) + min(nn))/2.;
        y_c = (max(mm) + min(mm))/2.;
        [pp,qq] = size(pixel);
        [coord_x,coord_y] = meshgrid(1:qq,1:pp);
        pixel_x = E.*((double(coord_x) - x_c).*pixel_scale).^2.*pixel_A;
        pixel_y = E.*((double(coord_y) - y_c).*pixel_scale).^2.*pixel_A;
        EIx(FN) = sum(sum(pixel_x));
        EIy(FN) = sum(sum(pixel_y));
    end
    [r1,p1]=corr(CSMI_outplane,EIx')
    [r2,p2]=corr(CSMI_inplane,EIy')
    [r3,p3]=corr(CSV_CSA,EA')
    [r4,p4]=corr(CSV_CSA,GA')
elseif par == 2
    for i_case = 1:10 %filenumber
        CTXA = dicomread(filelist{i_case});
        BB = Edge_Detect(CTXA);
        [by,bx] = find(BB==1);
        %outline = sort_smooth_outline([bx by]);
        figure(i_case);
        imshow(CTXA); hold on; axis equal; axis off;
        plot(bx,by,'r.'); pause;
    end;
end;







