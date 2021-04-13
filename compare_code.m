%% A sample code to generate Fig 1 and Table 1 of https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9303765
% This can be run once all the data is stored in the mat files, after running main.m  
% Monimoy Bujarbaruah

%%
clear all; 
close all; 
clc;
rng(2)

%% Load, plot and compare 
% Loading the existing data. Replace this data with new data after running the main.m file until convergence

% first cvx hull based method
load('allDataCVXHull.mat'); 
[prob_failcvx, cost_averagecvx] = MC_sim(XhatcvxOut, 100); 

% svm with eps = 0.3
load('allDataSVM30.mat');
XhatsvmOut30 = XhatsvmOut;
[prob_failsvm30, cost_averagesvm30] = MC_sim(XhatsvmOut30, 100); 

% svm with eps = 0.5
load('allDataSVM50.mat'); 
XhatsvmOut50 = XhatsvmOut;
[prob_failsvm50, cost_averagesvm50] = MC_sim(XhatsvmOut50, 100); 

%% Computing the true cost (if all constraints were known) 
[prob_failXnew, cost_averageXnew] = MC_sim(Xnew, 100); 

%% Plotting for Fig 1 of https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9303765
figure; 
plot(XhatsvmOut50,'color','r','alpha',0.9); hold on; 
hold on; 
plot(XhatsvmOut30,'color','b','alpha',0.5); grid on;
hold on; 
plot(Xnew,'color','k','alpha',0, 'linewidth',4);
legend({'With SVM: $\epsilon = 0.5$',...
    'With SVM: $\epsilon=0.3$','True','CVX Hull'},'Interpreter','latex','fontsize',25);
set(gca,'Fontsize',25,'fontweight','Bold');
