function G=LawyersData()
%addpath(append(pwd,'\datasets\1_Lawyers_data\'))
load('Lawyers.mat')
G=graph(Lawyers);

%Plot
% figure(1);
% plot(G)
% title('Original graph - Lawyers');
end