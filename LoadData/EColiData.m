function G=EColiData()
load('EColi.mat')
G=graph(EColi);
figure(1);
% plot(G);
% title('Original graph - EColi');
end