function x_guess = correctMapV2(ef, xf)

x_guess = xf;
tempI = find(ef(:,1) == 0);
x_guess(tempI) = x_guess(tempI) + ef(tempI,2);