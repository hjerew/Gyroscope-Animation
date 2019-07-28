function R = rMat(axis, sign, input)
input2=str2sym(input);
sign=lower(sign);   
    switch lower(axis)
        case 'x'
            if (strcmpi('p',sign))
                R = [1 0 0; 0 cos(input2) -sin(input2); 0 sin(input2) cos(input2)];
            else
                R = [1 0 0; 0 cos(input2) sin(input2); 0 -sin(input2) cos(input2)];
            end
        case 'y'
            if (strcmpi('p',sign))
                R = [cos(input2) 0 sin(input2); 0 1 0; -sin(input2) 0 cos(input2)];
            else
                R = [cos(input2) 0 -sin(input2); 0 1 0; sin(input2) 0 cos(input2)];
            end
        case 'z'
            if (strcmpi('p',sign))
                R = [cos(input2) -sin(input2) 0; sin(input2) cos(input2) 0; 0 0 1];
            else 
                R = [cos(input2) sin(input2) 0; -sin(input2) cos(input2) 0; 0 0 1];
            end
    end
end