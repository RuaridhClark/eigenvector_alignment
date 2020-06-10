% Counts the number of connected triples in a graph
% INPUTs: adjacency matrix
% OUTPUTs: integer - num conn triples
% Other routines used: kneighbors.m, loops3.m
% Note: works for undirected graphs only
% GB, Last updated: October 9, 2009
%
% Copyright (c) 2011, Massachusetts Institute of Technology. All rights
% reserved. Redistribution and use in source and binary forms, with or
% without modification, are permitted provided that the following
% conditions are met:
% 
% - Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer. - Redistributions
% in binary form must reproduce the above copyright notice, this list of
% conditions and the following disclaimer in the documentation and/or other
% materials provided with the distribution. - Neither the name of the
% Massachusetts Institute of Technology nor the names of its contributors
% may be used to endorse or promote products derived from this software
% without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function c=num_conn_triples(adj)

c=0;  % initialize

for i=1:length(adj)
    neigh=kneighbors(adj,i,1);
    if length(neigh)<2; continue; end  % handle leaves, no triple here
    c=c+nchoosek(length(neigh),2);
end

c=c-2*loops3(adj); % due to the symmetry triangles repeat 3 times in the nchoosek count