function G=PBData()
load('PB.mat')
G=graph(PB);
%Plot
% figure(1);
% plot(G);
% title('Original graph - Political Blogs');
end