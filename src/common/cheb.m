%/* ************************************************** */
%     Copyright 2014 Arash Bakhtiari
%
%     you may not use this file except in compliance with the License.
%     You obtain a copy of the License in the LICENSE file
%
%     Unless required by applicable law or agreed to in writing, software
%     distributed under the License is distributed on an "AS IS"" BASIS,
%     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%     See the License for the specific language governing permissions and
%     limitations under the License.
%/* ************************************************** */

classdef cheb < handle
%UNTITLED Summary of this class goes here
%   Detailed explanation goes here

    properties
    end

    methods (Static)

        %/* ************************************************** */
        function [x] = chebnodes1(N, xmin, xmax)
            if nargin < 2, xmin=-1; xmax=1; end;
            i = [0:N]';
            global CHEB_KIND
            if CHEB_KIND == 1
                x = -cos((i+1/2)*pi/(N+1));
            elseif CHEB_KIND == 2
                x = -cos(i*pi/N);
            end
            if nargin < 2, return; end;
            x = 0.5*(xmin+xmax) + 0.5*(xmax-xmin).*x;
        end

        %/* ************************************************** */
        function [xx, yy] = chebnodes2(resX, resY, xmin, xmax, ymin, ...
                                        ymax)
            if nargin < 2, resY = resX; end;
            if nargin < 3, xmin=-1; xmax=1; ymin=-1; ymax=1;end;
            x = cheb.chebnodes1(resX, xmin, xmax);
            y = cheb.chebnodes1(resY, ymin, ymax);
            [xx, yy] = meshgrid(x,y);
        end

        function [w] = chebcoeff(fn_val)
            global CHEB_KIND
            n1 = size(fn_val, 2);
            n2 = size(fn_val, 1);
            x = cheb.chebnodes1(n1-1);
            y = cheb.chebnodes1(n2-1);
            Tx = cheb.chebpoly(n1-1,x);
            Ty = cheb.chebpoly(n2-1,y);

            if CHEB_KIND == 1
                Nx = n1; 
                Ny = n2;
            elseif CHEB_KIND == 2
                Nx = n1-1;
                Ny = n2-1;
                fn_val(1,:) = fn_val(1,:)/2;
                fn_val(:,1) = fn_val(:,1)/2;
                fn_val(end,:) = fn_val(end,:)/2;
                fn_val(:,end) = fn_val(:,end)/2;
            end

            % using discrete orthogonality of chebyshev polynomials to compute the
            % coefficients for approximation using chebyshev
            % polynomial basis.
            for i=1:n1
                w_(:,i) = (reshape(fn_val(:,i),1,[])*Ty)*2/Ny;
            end
            for j=1:n2
                w(j,:) = (reshape(w_(j,:),1,[])*Tx)*2/Nx;
            end
        end

        %/* ************************************************** */
        function [fval] = chebeval2(w,x,y)
        % CHEBEVAL2(W, X, Y) Compute the values of a chebyshev
        % approximation at a regular grid specified by X, Y. where
        % W is the corresponding chebyshev coefficients.
            global CHEB_KIND
            fval = zeros(length(y),length(x));
            n1 = size(w,2);
            n2 = size(w,1);
            T_x = cheb.chebpoly(n1-1, x);
            T_y = cheb.chebpoly(n2-1, y);

            w(1,:) = w(1,:)/2;
            w(:,1) = w(:,1)/2;
            if CHEB_KIND == 2
                w(end,:) = w(end,:)/2;
                w(:,end) = w(:,end)/2;
            end

            for i=1:size(w,2)
                f_(:,i)=T_y*reshape(w(:,i),[],1);
            end
            for j=1:length(y)
                fval(j,:)=T_x*reshape(f_(j,:),[],1);
            end
        end

        %/* ************************************************** */
        function [T] = chebpoly(n,x)
        % CHEBPOLY(N, X) Retruns values of all chebyshev polynomials upto degree N
        % at points X
            x_ = reshape(x,[],1);
            T0 = ones(size(x_));
            T(:,1) = T0;
            if n == 0
                return;
            end
            T1 = x_;
            T(:,2) = T1;
            if n == 1
                return;
            end
            for i = 2:n
                T(:,i+1) = 2*x_.*T1-T0;
                T0 = T1;
                T1 = T(:,i+1);
            end
        end
    end
end
