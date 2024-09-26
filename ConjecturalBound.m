
LongRoots := function(type)
    roots := Roots(RootSystem(ReflectionGroup(type)));
    longroots := [];
    for i in [1..#roots] do 
        if IsLongRoot(RootSystem(ReflectionGroup(type)),i) eq true then 
            longroots := longroots cat [Transpose(Matrix(roots[i]))];
        end if;
    end for;
    return longroots;
end function;

InnerProd := function(x,y,Ematrix)
    m := #Rows(Ematrix);
    x := [x[i]/1 : i in [1..#x]];
    y := [y[i]/1 : i in [1..#y]];
    x := Ematrix*(Transpose(Matrix([x])));
    y := Ematrix*(Transpose(Matrix([y])));
    return Eltseq(&+[x[i]*y[i] : i in [1..m]])[1];
end function;

coroot := function(x,Ematrix);
    i := InnerProd(x,x,Ematrix);
    return [2*x[j]/i: j in [1..#x]];
end function;

PhiJ := function(posroots,J,n)
    nJ := {1..n} diff J;
    if #nJ eq 0 then 
        return posroots;
    else
        phij := [x : x in posroots | {x[i] : i in nJ} eq {0}];
        return phij;
    end if;
end function;

rhoJ := function(posroots,J,longroots,n)
    long := [Eltseq(x) : x in longroots];
    short := [x : x in posroots | not x in long];
    phij := PhiJ(posroots,J,n);
    shortj := [x : x in short | x in phij];
    longj := [x : x in long | x in phij];
    if #shortj eq 0 then 
        rhoJdash := [0 : i in [1..n]];
    else 
        rhoJdash := [&+[x[i] : x in shortj] : i in [1..n]];
    end if;
    if #longj eq 0 then 
        rhoJ := [0 : i in [1..n]];
    else
        rhoJ := [&+[x[i] : x in longj] : i in [1..n]];
    end if;
    return rhoJ,rhoJdash;
end function;

printinner := function(posroots,J,Ematrix,longroots,n)
    p,pdash := rhoJ(posroots,J,longroots,n);
    set := [];
    for x in posroots do 
        set := set cat [<x,InnerProd(coroot(x,Ematrix),pdash,Ematrix),InnerProd(coroot(x,Ematrix),p,Ematrix)>];
    end for;
    return set;
end function;

valpha := function(posroots,sign,J,Ematrix,longroots,n)
    innerset := printinner(posroots,J,Ematrix,longroots,n);
    valphaset := [[sign[1]*Numerator(x[2]),sign[2]*Numerator(x[3])] : x in innerset];
    return valphaset;
end function;

posorneg := function(coeff)
    if coeff eq [0,0] then
        return <[1,1],0>;
    elif (coeff[1] lt 0 and coeff[2] le 0) or (coeff[1] le 0 and coeff[2] lt 0) then 
        return <[-1,-1],0>;
    elif (coeff[1] gt 0 and coeff[2] ge 0) or (coeff[1] ge 0 and coeff[2] gt 0) then 
        return <[1,1],0>;
    elif coeff[1] lt 0 and coeff[2] gt 0 then 
        r := -coeff[2]/coeff[1];
        return <[1,-1],r>;
    elif coeff[1] gt 0 and coeff[2] lt 0 then 
        r := -coeff[2]/coeff[1];
        return <[-1,1],r>;
    else
        return "ERROR";
    end if ;
end function;

sum := function(x,y);
    return [x[1]+y[1],x[2]+y[2]];
end function;

nextstepposneg := function(coeff,weight)
    negweight := [-weight[1],-weight[2]];
    firststep := posorneg(coeff);
    negcoeff := [-coeff[1],-coeff[2]];
    nextstep := [];
    for i in [1..2] do 
        if firststep[1][i] eq 1 then 
            next := posorneg(sum(coeff,negweight));
            deg := [];
            for j in [1..2] do 
                if next[1][j] eq 1 then 
                    deg := deg cat [weight];
                elif next[1][j] eq -1 then 
                    deg := deg cat [coeff];
                else
                    return "ERROR";
                end if;
            end for;
            nextstep := nextstep cat [<deg,next[2]>];
        elif firststep[1][i] eq -1 then 
            next := posorneg(sum(negcoeff,negweight));
            deg := [];
            for j in [1..2] do 
                if next[1][j] eq 1 then 
                    deg := deg cat [weight];
                elif next[1][j] eq -1 then 
                    deg := deg cat [negcoeff];
                else
                    return "ERROR";
                end if;
            end for;
            nextstep := nextstep cat [<deg,next[2]>];
        else
            return "ERROR";
        end if;
    end for;
    return <firststep, nextstep>;
end function;

allranges := function(posroots,sign,J,Ematrix,longroots,n)
    long := [Eltseq(x) : x in longroots];
    short := [x : x in posroots | not x in long];
    steps := [];
    options := valpha(posroots,sign,J,Ematrix,longroots,n);
    for i in [1..#posroots] do 
        if posroots[i] in long then 
            weight := [0,2];
        elif posroots[i] in short then 
            weight := [2,0];
        else    
            return "ERROR";
        end if;
        steps := steps cat [nextstepposneg(options[i],weight)];
    end for;
    return steps;
end function;

choose := function(x,range);
    l := range[1];
    k := range[2];
    if l ge x[1][2] then 
        step := x[2][2];
    elif k le x[1][2] then 
        step := x[2][1];
    else
        "ERROR";
    end if;
    if l ge step[2] then 
        deg := step[1][2];
    elif k le step[2] then 
        deg := step[1][1];
    else
        "ERROR";
    end if;
    return deg;
end function;

degree := function(deg,rep)
    return deg[1]*rep[1]+deg[2]*rep[2];
end function;

//Doesn't work for F4 but does for G2 (figure out)
sumandbounds := function(posroots,rep,sign,J,Ematrix,longroots,lwo,n)
    all := allranges(posroots,sign,J,Ematrix,longroots,n);
    bounds := [y/1 : y in {x[1][2] : x in all} join {x[2][1][2] : x in all} join {x[2][2][2] : x in all}] cat [100];
    Sort(~bounds, func<u,v | u - v>);
    ranges :=[];
    for i in [1..(#bounds-1)] do 
        ranges := ranges cat [[bounds[i], bounds[i+1]]];
    end for;
    rangesdegrees := [];
    for x in ranges do 
        new := [];
        for y in all do 
            new := new cat [degree(choose(y,x),rep)];
        end for;
        rangesdegrees := rangesdegrees cat [lwo - (1/2)*(&+[x : x in new])];
    end for;
    return ranges, rangesdegrees;
end function;

conjecturalbound := function(type,J,sign)
    n := StringToInteger(type[2..#type]);
    root_system := RootSystem(ReflectionGroup(type));
    posroots := [Eltseq(Matrix(r)) : r in PositiveRoots(root_system)];
    longroots :=LongRoots(type);
    Wfin := CoxeterGroup(GrpFPCox, type);
    Q := Rationals();
    LL := PolynomialRing(Q);
    AssignNames(~LL, ["a"]);
    a := LL.1;
    L := PolynomialRing(LL);
    AssignNames(~L, ["b"]);
    b := L.1;
    simpleroots := [posroots[i] : i in [1..n]];
    aindices := [Index(posroots,x) : x in simpleroots | not x in longroots];
    bindices := [Index(posroots,x) : x in simpleroots | x in longroots];
    lwo := Eltseq(LongestElement(Wfin));
    if type eq "F4" then 
        Ematrix := Matrix([[0,0,0,1/2],[1,0,0,-1/2],[-1,1,0,-1/2],[0,-1,1,-1/2]]);
        rep := [a,b];
    elif type eq "G2" then 
        Ematrix := Matrix([[1/1,-2/1],[-1/1,1/1],[0,1/1]]);
        rep := [a,b];
    elif type[1] eq "A" then 
        row := [];
        for i in [1..n] do 
            row := row cat [[0/1 : j in [1..i-1]] cat [1/1,-1/1] cat [0/1 : j in [1..(n+1-i-1)]]];
        end for;
        Ematrix := Transpose(Matrix(row));
        rep := [1,1];
    elif type[1] eq "B" then 
        row := [];
        for i in [1..n-1] do 
            row := row cat [[0/1 : j in [1..i-1]] cat [1/1,-1/1] cat [0/1 : j in [1..(n-i-1)]]];
        end for; 
        row := row cat [[0/1 : j in [1..n-1]] cat [1]];
        Ematrix := Transpose(Matrix(row));
        rep := [a,b];
    elif type[1] eq "C" then 
        row := [];
        for i in [1..n-1] do 
            row := row cat [[0/1 : j in [1..i-1]] cat [1/1,-1/1] cat [0/1 : j in [1..(n-i-1)]]];
        end for; 
        row := row cat [[0/1 : j in [1..n-1]] cat [2]];
        Ematrix := Transpose(Matrix(row));
        rep := [a,b];
    elif type[1] eq "D" then 
        row := [];
        for i in [1..n-1] do 
            row := row cat [[0/1 : j in [1..i-1]] cat [1/1,-1/1] cat [0/1 : j in [1..(n-i-1)]]];
        end for; 
        row := row cat [[0/1 : j in [1..n-2]] cat [1,1]];
        Ematrix := Transpose(Matrix(row));
        rep := [1,1];
    elif type eq "E8" then 
        Ematrix := Transpose(Matrix([[1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2],[1,1,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],[0,0,-1,1,0,0,0,0],[0,0,0,-1,1,0,0,0],[0,0,0,0,-1,1,0,0],[0,0,0,0,0,-1,1,0]]));
        rep := [1,1];
    elif type eq "E7" then 
        Ematrix := Transpose(Matrix([[1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2],[1,1,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],[0,0,-1,1,0,0,0,0],[0,0,0,-1,1,0,0,0],[0,0,0,0,-1,1,0,0]]));
        rep := [1,1];
    elif type eq "E6" then 
        Ematrix := Transpose(Matrix([[1/2,-1/2,-1/2,-1/2,-1/2,-1/2,-1/2,1/2],[1,1,0,0,0,0,0,0],[-1,1,0,0,0,0,0,0],[0,-1,1,0,0,0,0,0],[0,0,-1,1,0,0,0,0],[0,0,0,-1,1,0,0,0]]));
        rep := [1,1];
    else
        return "Enter a type An,Bn,Cn,Dn,F4,G2,E6,E7 or E8";
    end if;
     lwo := rep[1]*#[x : x in lwo | x in aindices] + rep[2]*#[x : x in lwo | x in bindices];
    return sumandbounds(posroots,rep,sign,J,Ematrix,longroots,lwo,n);
end function;


    