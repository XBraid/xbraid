function cycleplot(filename);

X = load(filename);
Y = X(:,1)-X(:,2);
plot(Y);
view([0,-90]);
