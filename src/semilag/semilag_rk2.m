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

function [ cnew ] = semilag_rk2(xx,yy,zz,fconc_interp,fvel_interp,t)
%ADVECT_SL_RK2 Advect the c values one time step by using semi-lagrangian scheme
%
VPREVTSTEP  = 1;
VCURTSTEP   = 2;
VNEXTSTEP   = 3;
V2NEXTSTEP  = 4;

ti = t(VNEXTSTEP);
tf = t(VCURTSTEP);
n  = 10;

[xt,yt,zt] = trajectory_rk2(xx,yy,zz,fvel_interp,ti,tf,n);
cnew = fconc_interp(t(VCURTSTEP),xt,yt,zt);
end
