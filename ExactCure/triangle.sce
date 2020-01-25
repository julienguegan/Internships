/* triangle equilateral centre (x0,y0) et longueur l */

function triangle(x,l)
x0 = x(1)
y0 = x(2)

uy = y0 - sqrt(3)/6 ; ux = x0 - l/2
vy = y0 - sqrt(3)/6 ; vx = x0 + l/2
wx = x0 ; wy = y0 + sqrt(3)/3

plot(x0,y0,'k.')
plot([ux vx],[uy vy])
plot([ux wx],[uy wy])
plot([vx wx],[vy wy])

endfunction
