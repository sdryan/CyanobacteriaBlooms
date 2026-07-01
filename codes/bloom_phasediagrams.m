


A = load('phase_time_form.txt');
B = load('phase_time_decay.txt');
C = load('phase_time_overtake.txt');


size(A)
Mat_in = zeros(6,6);


for i = 1:size(A,1)
    if (A(i,4) > 119 && A(i,2) == .2)
        % x-comp f0, y-comp beta
        x1 = mod(A(i,3),2)+1;
        y1 = mod(A(i,1),.0025)+1;
        Mat_in(x1,y1) = 0;
    elseif (C(i,4) > 119 && A(i,2) == .2)
        % x-comp f0, y-comp beta
        x1 = mod(C(i,3),2)+1;
        y1 = mod(C(i,1),.0025)+1;
        Mat_in(x1,y1) = 1;
    else 
        % x-comp f0, y-comp beta
        x1 = mod(C(i,3),2)+1;
        y1 = mod(C(i,1),.0025)+1;
        Mat_in(x1,y1) = 2;
    end
end

x = linspace(0,10,6);
y = linspace(0,.0125,6);
[X,Y] = meshgrid(x,y);
[C,h] = contourf(X,Y,Mat_in,10)
clabel(C,h)