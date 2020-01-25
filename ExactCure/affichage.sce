clear
clf()
exec('C:\Users\Julien Guégan\Desktop\PFE\algorithmes\fonctions test.sce',-1)

/////  AFFICHAGE /////
function z = affiche(f,xmax,xmin,ymax,ymin,option)
    a = gcf()
    xstep = (xmax-xmin)/80
    axe_x = xmin:xstep:xmax
    ystep = (ymax-ymin)/80
    axe_y = ymin:ystep:ymax
    z=0
    if option=='2d' then
        for i = xmin:xstep:xmax
            ind_i = (i/xstep)-(xmin/xstep)+1
            z(ind_i) = f(u)
        end
        plot(x,z)
    elseif option=='3d'|option=='contour' then
        for i = xmin:xstep:xmax
            for j = ymin:ystep:ymax
                u = [i j]
                ind_i = (i/xstep)-(xmin/xstep)+1
                ind_j = (j/ystep)-(ymin/ystep)+1
                z(ind_i,ind_j) = f(u)
            end
        end
        if option == '3d' then
            plot3d1(axe_x,axe_y,z)
        elseif option == 'contour' then 
            xset("fpf"," ")//pour ne pas afficher les valeurs des courbes de niveau
            contour2d(axe_x,axe_y,z,50)
        end
    end
    a.color_map = rainbowcolormap(64)
endfunction
