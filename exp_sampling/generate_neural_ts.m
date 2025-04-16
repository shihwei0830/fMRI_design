function ts = generate_neural_ts(t_3Col,t_grid)
%
%
%
resel = t_grid(2)-t_grid(1);
nt = length(t_grid);
ts = zeros(nt,1);
n_ev = size(t_3Col,1);
for i=1:n_ev
    onset = t_3Col(i,1);
    [x,id] = find(t_grid == onset);
    dur_sec = t_3Col(i,2);
    dur_resel = dur_sec/resel;
    height = t_3Col(i,3);
    ts(id:id+dur_resel-1) = height;
end
