function G=CEleData()
load('CEle.mat')
G=graph(CEle);

%Plot
% figure(1);
% plot(G);
% title('Original graph - C.Elegans');
end