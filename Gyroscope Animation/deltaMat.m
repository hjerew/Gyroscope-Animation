function delta = deltaMat(r)
    delta=[r(2)^2 + r(3)^2 -r(1)*r(2) -r(1)*r(3); ...
        -r(2)*r(1) r(1)^2 + r(3)^2 -r(2)*r(3); ...
         -r(3)*r(1) -r(3)*r(2) r(1)^2 + r(2)^2];

end