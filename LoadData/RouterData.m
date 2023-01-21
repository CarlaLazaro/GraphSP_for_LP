function G=RouterData()
load('Router.mat')
G=graph(Router);
%Plot
% figure(1);
% plot(G);
% title('Original graph - Router');

end