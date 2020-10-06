figure;
% Load dataset
load('maneuvering_robot.mat');
truth = [TrueTrack.Trajectory.Vector];
obs = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',100^2,'Mapping',[1 3]);
%measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);

%% Cartesian
% subplot(1,2,1);
hold on;
truth_100 = truth.*100;
%h1 = plot(ax(1),measurement(1,k),measurement(2,k),'k*','MarkerSize', 10);
measurement = obs.feval(truth_100, true);
meas = obs.finv([measurement(1,:);measurement(2,:)]);
h1 = plot(meas(obs.Mapping(1),:),meas(obs.Mapping(2),:),'r+','MarkerSize', 15);
h2 = plot(truth_100(1,1:end),truth_100(3,1:end),'k--','LineWidth',2);
h3 = plot(truth_100(1,1),truth_100(3,1),'ko','MarkerSize', 20, 'MarkerFaceColor','Green');
h4 = plot(truth_100(1,end),truth_100(3,end),'ko','MarkerSize', 20, 'MarkerFaceColor','Red');
h5 = plot(0,0,'rd','MarkerSize', 20, 'MarkerFaceColor','Blue');
legend([h1,h2,h3,h4,h5],'Measurements','Trajectory', 'Start', 'End','Radar','Location','southeast')
str = sprintf('Vessel trajectory');
title(str)
xlabel('X (m)')
ylabel('Y (m)')
axis([0,2500,0,1500])
box on

%% Polar
% subplot(1,2,2);
% truth2meas = obs.heval(truth_100(1:2,:));
% % [a,b]=obs.heval(truth_100(1:2,:));
% %h1 = polarplot(a(x_start),b(y_start),'ko','LineWidth',2,'MarkerSize',15,'MarkerFaceColor','w');hold on;
% polarplot(truth2meas(1,:),truth2meas(2,:),'k-','LineWidth',2,'MarkerSize',10,'MarkerFaceColor','w');hold on;
% thetalim([0 90])
% %h2 = polarplot(a(x_end),b(y_end),'k^','LineWidth',2,'MarkerSize',15,'MarkerFaceColor','w');hold on;
% h3 = polarplot(measurement(1,:),100*measurement(2,:),'+r','MarkerSize',15);
% pax = gca;
% pax.ThetaAxisUnits = 'radians';
% %rlabel('Range (m)');
% %thetalabel('Bearing (rad)');