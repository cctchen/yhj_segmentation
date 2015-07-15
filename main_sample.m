clc; clear all; close all;

DXA = imread('sample.JPG');
DXA = rgb2gray(DXA);
bd_nodes = identify_boundary_pt_v2(DXA);
figure(1);
plot(bd_nodes(:,1),bd_nodes(:,2),'--'); axis equal; grid on;
   
