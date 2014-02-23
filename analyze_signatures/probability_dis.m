function probability_dis(data1, data2)
% Displays the probability distributions for each data set.

figure()
x_range =0:max(max(data1), max(data2));
sig1_N = hist(data1, x_range);
bar(x_range,sig1_N/trapz(x_range,sig1_N))

% Make see-through, prettier
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.75)
hold on

sig2_N = hist(data2, x_range);
bar(x_range,sig2_N/trapz(x_range,sig2_N))

% Make see-through
h = findobj(gca,'Type','patch');
set(h,'facealpha',0.75);

legend('data1', 'data2')
end