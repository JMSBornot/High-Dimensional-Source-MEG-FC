function [cluster, cluster_ext] = compute_cluster(E, ivert2irow, ivert2icol, ivertNeighBeg, indsel)

if (min(min(E)) < 1)
    error('<E> indices should start at 1, following Matlab notation.');
end
E = int32(E)-1;
if (min(ivert2irow) < 1)
    error('<ivert2irow> indices should start at 1, following Matlab notation.');
end
ivert2irow = int32(ivert2irow)-1;
if (min(ivert2icol) < 1)
    error('<ivert2icol> indices should start at 1, following Matlab notation.');
end
ivert2icol = int32(ivert2icol)-1;
if (min(ivertNeighBeg) < 1)
    warning('<ivertNeighBeg> indices should start at 1, following Matlab notation, at least that there are some disconnected points.');
end
ivertNeighBeg = int32(ivertNeighBeg)-1; % now indices in ivertNeighBeg start in 0 following C notation
if (min(indsel) < 1)
    error('<indsel> indices should start at 1, following Matlab notation.');
end
indsel = uint64(indsel)-1;
[cluster, cluster_ext] = compute_cluster_cppopt(E, ivert2irow, ivert2icol, ivertNeighBeg, indsel);
cluster = cluster+1;
nc = max(cluster); % number of clusters
cluster_ext = cluster_ext(1:nc);