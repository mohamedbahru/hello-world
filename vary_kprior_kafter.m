%% This script draws analytically the statistics as a function of kprior and kafter
function vary_kprior_kafter()
%% Setting the parameters of the model
model.dil_rate = log(2)/(3600*2); % dilution rate
model.krbs = 0.2; % Grieve 2005, Voigt1994 0.08 from CP
model.krnadecay = 0.0023 + model.dil_rate; % https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=111927
model.kd_prot = 1.9254e-04/1 + model.dil_rate; % https://bionumbers.hms.harvard.edu/bionumber.aspx?&id=108404
model.R = 1; % RNAp
model.CV2_param = 104; % constant that relates CV2 and mean
model.Sk_param = [-0.56, 18.37]; % pramaters that relates Skewness and mean

model.kbump = 1; % Inferred
model.b = 0.06; % Inferred

kafter = logspace(-4,1,100);
kprior = logspace(-4,-1,100);

[model.kafter, model.kprior] = meshgrid(kafter, kprior);

x = [30, 200]; % TSS distance 

%% Calculating Pmean
for i = 1:numel(x)
Out(i)= get_stats(model, x(i));
end

%% Making the surface plot

f1 = surface_plot(model.kafter, model.kprior, Out, 'PMean', ...
    [[0 max(model.kafter(:))]], [0 max(model.kprior(:))], [0 inf],...
    {'kafter', 'kprior', 'PMean'}, {'x=30', 'x=200'});

f2 = surface_plot(model.kafter, model.kprior, Out, 'PCV2', ...
    [[0 max(model.kafter(:))]], [0 max(model.kprior(:))], [0 inf],...
    {'kafter', 'kprior', 'PCV2'}, {'x=30', 'x=200'});

f3 = surface_plot(model.kafter, model.kprior, Out, 'PSkew', ...
    [[0 max(model.kafter(:))]], [0 max(model.kprior(:))], [0 inf],...
    {'kafter', 'kprior', 'PSkew'}, {'x=30', 'x=200'});

save_figure(f1, './figures', 'Pmean')
save_figure(f2, './figures', 'PCV2')
save_figure(f3, './figures', 'PSkew')

end
%% Function to save the figures. This will create and save the plots
function save_figure(f, pathname, filename)

% Creating the path in the results if not exists
if ~exist(pathname)
    mkdir(pathname)
end

savefig(f, fullfile(pathname,sprintf('%s.fig',filename)), 'compact')
saveas(f, fullfile(pathname,sprintf('%s.tif',filename)))

end
%% function to get the statics for certain x
function Out = get_stats(model, x)
dil_rate = model.dil_rate;
krbs = model.krbs;
krnadecay = model.krnadecay;
kd_prot = model.kd_prot;
R = model.R;
kbump = model.kbump;
b = model.b;
kafter = model.kafter;
kprior = model.kprior;
CV2_param = model.CV2_param;
Sk_param = model.Sk_param;

% Mean protein
PMean = 1./(kprior*R) + 1./kafter + kbump./(kprior*R.*kafter) * exp(-b*x) * 1./(1+(kafter./(kprior*R)));
PMean = 1./PMean * 2*krbs/(krnadecay*kd_prot);

% CV2 protein
PCV2 = CV2_param./PMean;

% Skewness protein
PSkew = Sk_param(2) * PMean.^Sk_param(1);

% Creating the output handle
Out.PMean = PMean;
Out.PCV2 = PCV2;
Out.PSkew = PSkew;

end

%% Function to make a surface plot
function fig = surface_plot(x_data, y_data, z_data, stat, x_bounds, y_bounds,...
z_bounds, labels, legends)

% x_data - data for x
% y_data - data for y
% z_data - data for z
% x_bounds - boundaries for x
% y_bounds - boundaries for y
% z_bounds - boundaries for z
% labels - {'x_label', 'y_label', 'z_label'}
% legends - legends for the surface plot

eval(sprintf('z = {z_data.%s};', stat));

fig = figure;

surf(x_data,y_data,z{1}, 'EdgeColor', 'none', 'FaceAlpha',0.5, 'FaceColor',[0.4 0 0]);
hold on;
surf(x_data,y_data,z{2}, 'EdgeColor', 'none', 'FaceAlpha',0.5, 'FaceColor',[0 0.4 0]);

% Set the remaining axes properties
set(gca,'XMinorTick','on','XScale','log','YMinorTick','on','YScale','log',...
    'ZMinorTick','on','ZScale','log');

xlabel(labels{1})
ylabel(labels{2})
zlabel(labels{3})

legend1 = legend(legends);
xlim(x_bounds)
ylim(y_bounds)
zlim(z_bounds)

% Create legend
set(legend1,...
    'Position',[0.729324406399434 0.88571428798494 0.151785712476288 0.0869047596341088]);
set(gca,'FontSize',13);
end