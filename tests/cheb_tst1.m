function cheb_tst1()
    clear all;
    close all;
    addpath('../src/common/');

    global CHEB_KIND;
    CHEB_KIND          = 1;

    f = @(x,y) sin(x).*cos(y);
    %f = @get_gaussian;

    resPerNode = 3;
    Nx = resPerNode;
    Ny = Nx;
    xg = cheb.chebnodes1(Nx);
    yg = cheb.chebnodes1(Ny);
    [xxg, yyg] = cheb.chebnodes2(Nx, Ny)

    fgval = f(xxg,yyg);
    w = cheb.chebcoeff(fgval);

    [xxg_fun, yyg_fun] = chebpts2(Nx+1, Ny+1);
    fgval_fun = f(xxg_fun,yyg_fun);
    ffun = chebfun2(fgval_fun);
    wfun = chebcoeffs2(ffun);

    fffun = chebfun2(f)

    xq = linspace(0,0.5,3);
    yq = xq;
    [xxq, yyq] = meshgrid(xq, yq);

    xq = [0.5 0.6];
    yq = [0.7 0.8];
    [xxq, yyq] = meshgrid(xq, yq);

    % xq = [0.5];
    % yq = [0.7];
    % xxq = [0.5];
    % yyq = [0.7];

    % xq = [0.99];
    % yq = [0.99];
    % xxq = [0.99];
    % yyq = [0.99];

    fe = f(xxq,yyq)
    ff = ffun(xxq, yyq)
    fff = fffun(xxq, yyq)
    fi = cheb.chebeval2(w, xq, yq)

    figure
    plot(ffun)
    figure
    plot(fffun)
    figure
    x = linspace(-1,1,100);
    y = x;
    [X,Y] = meshgrid(x, y);
    surf(X, Y, cheb.chebeval2(w, x, y))

    dfun = fe-ff;
    fprintf('CHEBFUN ERROR: ');
    err = max(abs(dfun(:)))

    dffun = fe-fff;
    fprintf('CHEBFUN FUNCTION ERROR: ');
    err = max(abs(dffun(:)))

    d = fe-fi;
    fprintf('MY ERROR: ');
    err = max(abs(d(:)))

    %/* ************************************************** */
    function value = get_gaussian(x,y)
        OM = 1;
        t = 0;
        xc = 0.5;
        yc = 0.5;
        xci = 0.6;
        yci = 0.5;

        [alpha,RHO] = cart2pol(xci-xc,yci-xc);
        alphat = alpha + t*OM;

        [xct,yct] = pol2cart(alphat,RHO);
        xct = xct + xc;
        yct = yct + yc;

        theta = 0;
        sigmax = 0.06;
        sigmay = 0.06;
        value = gaussian(x,y,xct,yct,theta,sigmax,sigmay);
    end
end