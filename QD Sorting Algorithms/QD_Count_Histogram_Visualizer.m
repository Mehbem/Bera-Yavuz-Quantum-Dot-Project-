% Bin Labels/Edges
Bins = 400:100:1000;

% Bin Counts
Bin_counts = [82 77 71 63 58 50 42];


% Plot histogram-like bar chart
figure()
bar(Bins, Bin_counts);

xlabel('Photon Count');
ylabel('Number of QD in Range');
for i = 1:numel(Bin_counts)
    label_text = sprintf("%d",Bin_counts(i));
    text(Bins(i),Bin_counts(i)+3,label_text,"FontWeight","bold")
end
title('QDs in the range 894.5-894.7 nm organized by arb. counts');