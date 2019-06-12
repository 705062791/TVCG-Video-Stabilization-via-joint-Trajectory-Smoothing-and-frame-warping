%% Show the feature point connection between two Images.
function result = shFtr(I1, I2, I1f, I2f)
imshow(I1);
height = max(size(I1, 1), size(I2, 1));
width = size(I1, 2);
channel = size(I1, 3);
n_ftr = size(I1f, 1);
if n_ftr ~= size(I2f, 1)
	printf('not match');
	exit;
end
cmap = hsv(6);
hold on
plot(I1f(:,1),I1f(:,2),'r.');
for ftrIndex = 1:n_ftr
	x1 = I1f(ftrIndex, 1);
	y1 = I1f(ftrIndex, 2);
	x2 = I2f(ftrIndex, 1);
	y2 = I2f(ftrIndex, 2);
	line([x1,x2],[y1,y2],'Color',cmap(mod(ftrIndex, 6)+1,:),'LineWidth',1);
	%if mod(ftrIndex, 10) == 0
	%	pause
	%	imshow([I1 I2]);
	
	%end
end