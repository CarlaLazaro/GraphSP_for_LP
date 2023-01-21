function G=YeastData()
load('Yeast.mat')
G=graph(Yeast);
% figure(1)
% plot(G);
% title('Original graph - YeastData');
end