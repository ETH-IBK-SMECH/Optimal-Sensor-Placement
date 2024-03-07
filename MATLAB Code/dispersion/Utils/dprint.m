function s = dprint(v,str,form)

%
% DPRINT  
%
% prints a column vector V with 
% a centered caption STR on top.
% 
% s=dprint(v,str,form)
%
% FORM argument is a string containing
% C language conversion specifications.
% See the help of SPRINTF command.
%
% See also FPRINTF, SPRINTF.
%
% Reference page in Help browser:
% <a href="matlab: web([docroot '/toolbox/mdac/funref/dprint.html'],'-helpbrowser')">doc dprint</a>
%

%
% Author:
% 1st Ed:
% Last Update:
% The Mathworks Corporation
% Copyright Mathworks
%

nv = size(v,1);
lstr = length(str);
lrpad = char(' '*ones(nv+4,2));  % left and right spacing
if isempty(v), 
   s = [];
   return
end
% Convert V to string (if not)
if isa(v,'char') == 1
    s = v;
else
    rev = real(v);
    s = [blanks(nv)' num2str(abs(rev),form)];
    s(rev<0,1) = '-';
    if ~isreal(v),
       % Add complex part
       imv = imag(v);
       imags = num2str(abs(imv),[form 'i']);
       imags(~imv,:) = ' ';
       signs = char(' '*ones(nv,3));
       signs(imv>0,2) = '+';
       signs(imv<0,2) = '-';
       s = [s signs imags];
    end
end
% Dimensions
ls = size(s,2);
lmax = max(ls,lstr);
ldiff = lstr - ls;
ldiff2 = floor(ldiff/2);
str = [blanks(-ldiff2) str blanks(-ldiff+ldiff2)];
s = [char(' '*ones(nv,ldiff2)) s char(' '*ones(nv,ldiff-ldiff2))];
% Put pieces together
s = [blanks(lmax) ; str ; blanks(lmax) ; s ; blanks(lmax)];
s = [lrpad s lrpad];













