function hcc_trace_matrix = Markov_Prob(TransArray)
trace_matrix = NaN(52,10);
trace_matrix(1,:) = [1,0,0,0,0,0,0,0,0,0];
trace_matrix(2,:) = [1,0,0,0,0,0,0,0,0,0]*TransArray;
for i = 3:52
    trace_matrix(i,:) = trace_matrix(i-1,:) * TransArray;
end
%hcc = half cycle corrected
hcc_trace_matrix = NaN(52,10);
hcc_trace_matrix(1,:) = 0.5*[1,0,0,0,0,0,0,0,0,0] + 0.5*trace_matrix(2,:);
for i = 2:51
    hcc_trace_matrix(i,:) = 0.5*trace_matrix(i,:) + 0.5*trace_matrix(i+1,:);
end
hcc_trace_matrix(52,:) = trace_matrix(52,:);
end