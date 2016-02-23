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
        function [x] = chebnodes1(resX, xmin, xmax)
            if nargin < 2, xmin=-1; xmax=1; end;
            n1 = resX+1;
            i = [0:(n1-1)]';
            x = -cos((i+1/2)*pi/n1);
            if nargin < 2, return; end;
            x = 0.5*(xmin+xmax) + 0.5*(xmax-xmin).*x;

            % i = [0:(n1-1)]';
            % x=cos((i+1/2)*pi/n1);
        end

        %/* ************************************************** */
        function [xx, yy] = chebnodes2(resX, resY, xmin, xmax, ymin, ...
                                        ymax)
            if nargin < 2, resY = resX; end;
            if nargin < 3, xmin=-1; xmax=1; ymin=-1; ymax=1;end;
            x = cheb.chebnodes1(resX, xmin, xmax);
            y = cheb.chebnodes1(resY, ymin, ymax);
            [xx, yy] = meshgrid(x,y);

            % n1 = resX+1;
            % n2 = resY+1;
            % xx=ones(n1,n2);
            % yy=ones(n1,n2);
            % for i=1:n1
            %     xx(i,:)=xx(i,:)*cos((i-1/2)*pi/n1);
            % end
            % for j=1:n2
            %     yy(:,j)=yy(:,j)*cos((j-1/2)*pi/n2);
            % end
        end

        function [w] = chebcoeff(fn_val)
            n2 = size(fn_val, 1);
            n1 = size(fn_val, 2);
            x = cheb.chebnodes1(n1-1);
            y = cheb.chebnodes1(n2-1);

            Tx = cheb.chebpoly(n1-1,x);
            Ty = cheb.chebpoly(n2-1,y);

            % using discrete orthogonality of chebyshev polynomials to compute the
            % coefficients for approximation using chebyshev
            % polynomial basis.
            for i=1:n1
                w_(:,i) = (reshape(fn_val(:,i),1,[])*Ty)*2/n2;
            end
            for j=1:n2
                w(j,:) = (reshape(w_(j,:),1,[])*Tx)*2/n1;
            end
            w(1,:) = w(1,:)/2;
            w(:,1) = w(:,1)/2;
        end

        %/* ************************************************** */
        function [fval] = chebeval2(w,x,y)
        % CHEBEVAL2(W, X, Y) Compute the values of a chebyshev
        % approximation at a regular grid specified by X, Y. where
        % W is the corresponding chebyshev coefficients.
            fval = zeros(length(y),length(x));
            n1 = size(w,2);
            n2 = size(w,1);
            T_x = cheb.chebpoly(n1-1, x);
            T_y = cheb.chebpoly(n2-1, y);
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
