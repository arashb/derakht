function cheb_tst1()
    clear all;
    close all;
    addpath('../src/common/');

    resPerNode = 30;
    xg = cheb.chebnodes1(resPerNode);
    yg = xg;
    [xxg, yyg] = cheb.chebnodes2(resPerNode, resPerNode)

    f = @(x,y) sin(x).*cos(y);
    % fgval = f(yyg,xxg)
    fgval = f(xxg,yyg)

    w = cheb.chebcoeff(fgval)

    ffun = chebfun2(fgval);
    wfun = chebcoeffs2(ffun)

    xq = linspace(0,0.5,3);
    yq = xq;
    [xxq, yyq] = meshgrid(xq, yq);

    % xq = [0.5 0.6];
    % yq = [0.7 0.8];
    % [xxq, yyq] = meshgrid(xq, yq)

    xq = [0.5];
    yq = [0.7];
    xxq = [0.5];
    yyq = [0.7];

    fe = f(xxq,yyq)
    ff = ffun(xxq, yyq)

    fi = cheb.chebeval2(w, xq, yq)
    % for xindx =1:length(xq)
    %     for yindx=1:length(yq)
    %         fi(yindx, xindx) = cheb.chebeval2(w, xq(xindx), yq(yindx));
    %     end
    % end
    fi

    dfun = fe-ff;
    err = max(abs(dfun(:)))

    d = fe-fi;
    err = max(abs(d(:)))
end