function f = mini(x, WLength, WRes,H)
       r = sqrt(2*x(2));
       f =  sum(sum(abs(CreateDisplacedThermalHusimiPhaseAveraged( x(1), WLength, WRes, r,'PlotOption',false)-H)));
    end