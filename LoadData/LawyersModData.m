function G=LawyersModData()
%addpath(append(pwd,'\datasets\1_Lawyers_data\'))
load('Lawyers_mod.mat')
G=H;

%Plot
% figure(1);
% plot(G)
% title('Original graph - Lawyers');
end