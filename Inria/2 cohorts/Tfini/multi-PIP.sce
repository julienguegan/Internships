clear
clf()

for J = 1:3
    for I = 1:3
        perdants = []
        gagnants = []

        M1σ = 2950+I*250
        M2σ = 2950+J*250

        exec('C:\Users\Julien Guégan\Documents\Cours\MAM4\STAGE\2 cohorts\Tfini\PIP - Tfini.sce',-1)

        subplot(3,3,(I-1)*3+J)

        //Mσ =  3000:3800
        //plot(Mσ,Mσ,'k')
        plot(perdants,gagnants,'r.')
        plot(gagnants,perdants,'b.')
       // contour(Mσ,Mσ,s,[0.99,0.999,1,1.001,1.01])
    end
end
