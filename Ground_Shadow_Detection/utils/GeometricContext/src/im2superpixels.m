function imsegs = im2superpixels(im)

prefix = num2str(floor(rand(1)*10000000));
fn1 = ['./tmpim' '.ppm'];
fn2 = ['./tmpimsp' '.ppm'];
segcmd = 'C:\Users\hasanbalci\Documents\MATLAB\Ground_Shadow_Detection\utils\GeometricContext\src\segment\segment 0.8 100 100';
%'/home/dhoiem/cmu/tools/superpixelsFH/segment 0.8 100 100';%

imwrite(im, fn1);
system([segcmd ' ' fn1 ' ' fn2]);
imsegs = processSuperpixelImage(fn2);

%delete(fn1);
%delete(fn2);