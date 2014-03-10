%Generate plotting mesh for Square Coax example
function [sol2D] = MyPlotMesh(innerBC,outterBC, x_vals, x, y, x_part, y_part)

% Define the partition steps for the three X and Y intervals.
for i=1:3
    x_inter(i) = x(i+1) - x(i);
    y_inter(i) = y(i+1) - y(i);
end

% Determine the step size for each interval of x, store it in a.
a(1) = x_inter(1)/x_part;
a(3) = x_inter(3)/x_part;
a(2) = floor((2*x_inter(2))/(a(1)+a(3)));
max_x_elm = 2*x_part + a(2) + 1;
x_subs = a(2);
a(2) = x_inter(2)/a(2);

% Determine the step size for each interval of y, store it in b.
b(1) = y_inter(1)/y_part;
b(3) = y_inter(3)/y_part;
b(2) = floor((2*y_inter(2))/(b(1)+b(3)));
max_y_elm = 2*y_part + b(2) + 1;
y_subs = b(2);
b(2) = y_inter(2)/b(2);

sol2D = zeros(max_x_elm);
count = 1;
for i=1:(max_x_elm-1)
    if(i==1 || i==max_x_elm)
        continue;
    end
    
    for j=1:max_y_elm
        if(j==1 || j==max_y_elm)
           continue;
        elseif((i>x_part && i<=(max_x_elm-x_part)) && (j>y_part && j<=(max_y_elm-y_part)))
            sol2D(i,j) = innerBC;
            continue;
        end
        sol2D(i,j) = x_vals(count);
        count = count + 1;
    end
end
