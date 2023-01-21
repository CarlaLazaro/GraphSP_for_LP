function G=USAirData()
load('USAir.mat')
G=graph(USAir);
%Plot
% figure(1);
% plot(G)
% title('Original graph - USAir');
end