# code to do network analysis of AP 
# genes in the context of Reactome.bd
# networks.
# Under construction.

# loads data into the workspace as global vars 
# c_exprs and c_condition
get_colon_data = function(){
  library(hgu133plus2.db)
  load("apData_gpl570.rda")
  idx = pData(apData_gpl570)$Tissue == "colon"
  c_exprs <<- exprs(apData_gpl570)[,idx]
  c_condition <<- pData(apData_gpl570)[idx,]$Status
}

# Takes expression data in, indexed by probe_ids
# then converts them to Entrez ids, 
# collapsing the multiple probes that match to one
# entrez ids via the median
convert_probes_entrez = function(expr_data){
  library(hgu133plus2.db)
  
  probe_ids = rownames(expr_data)
  entrez_ids = hgu133plus2ENTREZID[probe_ids]
  entrez_ids = toTable(entrez_ids)
  
  library(plyr)
  
  dfed = as.data.frame(expr_data)
  dfed = dfed[entrez_ids$probe_id,]
  
  dfed$gene = as.numeric(entrez_ids$gene_id)
  
  z_median = ddply(dfed, 'gene', function(z) apply(z, 2, median))
  
  rownames(z_median) = z_median[,'gene']
  z_median$gene = NULL
  
  return(as.matrix(z_median))
}

# alternative, possibly faster method
# untested
convert_probes_entrez2 = function(probeSetExpr){
  probe_ids = rownames(probeSetExpr)
  entrez_ids = hgu133plus2ENTREZID[probe_ids]
  egId = toTable(entrez_ids)
  
  #probeSetExpr # matrix nprobes by nsamples
  #egId # vector nprobes with ngene unique ids
  egIndex = split(seq(along=egId), egId)
  geneExpr = sapply(egIndex, function(ind) matrixStats::colMedians(probeSetExpr[ind,,drop=FALSE]))
  
  return(geneExpr)
}

convert_ap_probes = function(ap){
  library(hgu133plus2.db)
  probes = getProbesetIds(ap)
  gids = toTable(hgu133plus2ENTREZID[getProbesetIds(ap)])
  gids = gids$gene_id
  gids = unique(gids)
  return(gids)
}

# alternate method, using aggregate instead of ddply
# under construction
# does not work
convert_probes_2 = function(expr_data){
  library(hgu133plus2.db)
  
  probe_ids = rownames(expr_data)
  entrez_ids = hgu133plus2ENTREZID[probe_ids]
  entrez_ids = toTable(entrez_ids)
  
  median_data = expr_data[entrez_ids$probe_id]
  #median_data = transform(median_data, eid=eids$gene_id)
  
  median_data = as.data.frame(median_data)
  show(median_data)
  median_data = aggregate(median_data, by=list(entrez_ids$gene_id), FUN=median, na.rm=FALSE)
  
  return(median_data)
}

# data exploration
gene_correlation = function(data, threshold=.9){
  ed = exprs(data)
  status = pData(data)$Status == 1
  
  control_samples = ed[,status]
  cancer_samples = ed[,!status]
  
  control_cors = cor(t(control_samples))
  control_cors = control_cors - diag(1, dim(control_cors)[1])
  control_cors = control_cors[lower.tri(control_cors)]
  
  cancer_cors = cor(t(cancer_samples))
  cancer_cors = cancer_cors - diag(1, dim(cancer_cors)[1])
  cancer_cors = cancer_cors[lower.tri(cancer_cors)]
  
  high_cors = control_cors
  high_cors = high_cors[control_cors > .9]
  
  high_ind = control_cors > threshold
  low_ind = control_cors < -threshold
  
  #plot(control_cors[high_ind], cancer_cors[high_ind])
  plot(control_cors[low_ind], cancer_cors[low_ind])
}


# builds anti-profile...
# useful functions:
# getProbesetIds(ap)
# getNormalRegions(ap)
build_ap = function(data, status){
  library(antiProfiles)
  stats = apStats(data, status)
  ap = buildAntiProfile(stats, tissueSpec=FALSE)
  return(ap)
}

get_entrez_ids = function(data){
  x = hgu133plus2ENTREZID[rownames(data)]
  return(toTable(x))
}

get_ensemble_ids = function(data){
  x = hgu133plus2ENSEMBL[rownames(data)]
  return(toTable(x))
}

write_entrez_ids = function(tab, f="eids.txt"){
  write(tab$gene_id, file=f)
}

write_ensembl_ids = function(tab, f="eids.txt"){
  write(tab$ensembl_id, file=f)
}

# loads relations from the file rels.txt
# which is created by the python code
load_gene_relations = function(data, rfile, f="rels.txt"){
  
  r = read.table(f)
  m = data
  m[,] = 0
  
  for (i in 1:nrow(r)){
    m[as.character(r[i,1]), as.character(r[i,2])] = 1
  }
  
  return(m)
}


# old
rel_from_pathid_broad = function(pathid, data){
  r = toTable(reactomePATHID2EXTID[pathid])$gene_id
  
  m = mat.or.vec(nrow(data), nrow(data))
  rownames(m) = rownames(data)
  colnames(m) = rownames(data)
  
  r = intersect(rownames(data), r)
  
  for (i in 1:length(r)){
    for (j in i:length(r)){
      m[which(rownames(m) == r[i]),which(rownames(m)==r[j])] = 1
    }
  }
  
  return(m)
}

# new
rel_from_pathid_broad2 = function(pathid, data){
  r = toTable(reactomePATHID2EXTID[pathid])$gene_id
  
  m = mat.or.vec(nrow(data)*(nrow(data)-1), 2)
  #rownames(m) = rownames(data)
  #colnames(m) = rownames(data)
  
  r = intersect(rownames(data), r)
  
  rn = rownames(data)
  
  c = 1
  
  for (i in 1:(length(rn)-1)){
    for (j in (i+1):length(rn)){
      if (!(rn[i] %in% r) & !(rn[j] %in% r)){
        m[c,1] = rn[i]
        m[c,2] = rn[j]
        c = c + 1
      }
    } 
  }
  
  return(m[1:c-1,])
}

# old version
rel_from_pathid_slim = function(pathid, data){
  r = toTable(reactomePATHID2EXTID[pathid])$gene_id
  r = intersect(rownames(data), r)
  
  m = mat.or.vec(length(r)*length(r),2)
  
  write(r, 'eids.txt')
  system("python extract_reactome.py eids.txt rels.txt")
  
  show(r)
  s_data = data[r,]
  s_data = var(t(s_data))
  rels <<- load_gene_relations(s_data)
  rels <<- rels + t(rels)
  c = 1;
  
  for (i in 2:length(r)){
    for (j in 1:(i-1)){
      show(c(i,j))
      if (rels[i,j] == 0){
        #m[c,] = c(which(rownames(m) == r[i]), which(rownames(m) == r[j]))
        m[c,] = c(i,j)
        c = c+1
      }
    }
  }
  
  return(m[1:c-1,])
}

# working version
# TODO: document
rel_from_pathid_slim2 = function(pathid, data){
  r = toTable(reactomePATHID2EXTID[pathid])$gene_id
  r = intersect(rownames(data), r)
  
  m = mat.or.vec(length(r)*length(r),2)
  
  write(r, 'eids.txt')
  system("python extract_reactome.py eids.txt rels.txt")
  
  #show(r)
  s_data = data[r,]
  s_data = var(t(s_data))
  rels = load_gene_relations(s_data)
  rels = rels + t(rels)
  c = 1;
  
  for (i in 2:length(r)){
    for (j in 1:(i-1)){
      #show(c(i,j))
      if (rels[i,j] == 0){
        #m[c,] = c(which(rownames(m) == r[i]), which(rownames(m) == r[j]))
        m[c,] = c(i,j)
        c = c+1
      }
    }
  }
    
  return(list(rels, m[1:c-1,]))
}


# gets the subset of the expression data
# for all genes within the given pathway
d_sub = function(pathid, data){
  r = toTable(reactomePATHID2EXTID[pathid])$gene_id
  r = intersect(rownames(data), r)  
  return(data[r,])
}


# plots the differences between the partial correlations
# for the cancer and normal samples
# controlg, cancerg = output from glasso
# rels = adjacency matrix for network
# ap = list of ap genes entrez ids
plot_differences = function(controlg, cancerg, rels, ap){
  i = which(rels != 0, arr.ind = T)
  #show(i)
  #show(controlg$wi)
  #show(cancerg$wi)
  x = controlg$wi[i]
  y = cancerg$wi[i]
  
  #show(rels)
  #show(ap)
  
  i = which(rels[ap,] != 0, arr.ind = T)
  x2 = controlg$wi[i]
  y2 = cancerg$wi[i]
  
  #show(ap)
  #show(i)
  #show(x2)
  
  plot(x, y, main="glasso Fitted Partial Correlation Coefficients", xlab="Control Samples", ylab="Cancer Samples", pch=1, col=1)
  #  plot(c(x, x2), c(y, y2), main="glasso Fitted Partial Correlation Coefficients", 
  #       xlab="Control Samples", ylab="Cancer Samples", 
  #       pch=1, col=1)
  points(x2, y2, col=2)
  
}

# loads the colon data
# prints gets several pathways, and runs glasso on them
# prints the results to a chart
sample_paths_glasso = function(){
  library(antiProfiles)
  library(glasso)
  library(hgu133plus2.db)
  library(reactome.db)
  
  #get_colon_data()	# should use setup_colon_ap_data instead

  # pathways list:
  # apoptosis, degredation of the extracellular matrix, growth hormone receptor signaling (immune -> cytokine),
  # signaling by EGFR, Signalling by NOTCH1
  sample_paths = c("169911", "1474228", "982772", "177929", "2644602")
  
  library(reactome.db)
    
  # get ap stuff...  # should use setup_colon_ap_data instead
  #ap = build_ap(c_exprs, c_condition)
  #ap_entrez_ids = convert_ap_probes(ap)
  
  control_samples = c_condition == 0
  cancer_samples = c_condition == 1
  
  
  #for (i in 1:length(sample_paths)){   # TODO: why doesn't this work?
  for (i in 5:5){
    data = load_path_exprs_subset(sample_paths[i])
    data = convert_probes_entrez(data)
    r <<- rel_from_pathid_slim2(sample_paths[i], data)
    rel <<- r[[1]]
    
    zero_partial_corrs <<- r[[2]]
    dapgl = ap_entrez_ids
    dapgl = dapgl %in% rownames(data)
    dapgl = ap_entrez_ids[dapgl]
    print("subset found")

    cov_control = data[,control_samples]
    cov_cancer = data[,cancer_samples]
    print("subsets made")

    cov_control = var(t(cov_control))
    cov_cancer = var(t(cov_cancer))
    print("cov done")
    
    
    controlg = glasso(cov_control, rho=0.0001, zero=zero_partial_corrs)    
    cancerg = glasso(cov_cancer, rho=0.0001, zero=zero_partial_corrs)
    show("glasso done")

    plot_differences(controlg, cancerg, rel, dapgl)
  }
}


# Finds and returns all pathways that contain genes from the
# ap_entrez_genes list
# under construction...
find_pathways = function(ap_entrez_genes){
  pathids = c()
  for (i in 1:length(ap_entrez_genes)){
    #  	if (ap_entrez_genes[i] %in% reactomeEXTID2PATHID){
    r = reactomeEXTID2PATHID[ap_entrez_genes[i]]
    length(r)
    r = toTable(r)
    r = r$DB_ID
    pathids = c(pathids, r)
    #    }
  }
  
  return(unique(pathids))
}

# returns the subset of the colon exprs object which has
# probe ids that map to genes in the pathway defined by
# pathid
load_path_exprs_subset = function(pathid){
  probe_ids = rownames(c_exprs)
  r = toTable(reactomePATHID2EXTID[pathid])$gene_id
  c = toTable(hgu133plus2ENTREZID[probe_ids])
  sub = c$probe_id[c$gene_id %in% r]
  return(c_exprs[sub,])
}

# loads all of the required data for ap network analysis
# gets colon_data, and it builds the anti-profile
# translates the ap probeset list to entrez ids
setup_colon_ap_data = function(){
  get_colon_data()
  ap <<- build_ap(c_exprs, c_condition)
  ap_entrez_ids <<- convert_ap_probes(ap)
}

# gets an adjacency matrix for a reactome.db network
# corresponding to a gene pathway, where all genes in
# the network must be part of the c_exprs object
get_graph_am = function(pathid){
  data = load_path_exprs_subset(pathid)
  data = convert_probes_entrez(data)
  r = rel_from_pathid_slim2(pathid, data)
  rel = r[[1]]
  
  return(rel)
}

# gets some basic AP/Non-AP statistics for the network  
# am = adjacency matrix, with entrez ids as the row and col names
# ap_genes = list of entrez ids
graph_ap_stats = function(am, ap_genes){
  library(igraph)
  g = graph.adjacency(am)
  idx = rownames(am) %in% ap_genes
  rs = rowSums(am)
  show("Degree of Non-AP Genes:")
  show(rs[!idx])
  show("Degree of AP Genes:")
  show(rs[idx])
  show("Avg Non-AP Gene Degree:")
  show(mean(rs[!idx]))
  show("Avg AP Gene Degree:")
  show(mean(rs[idx]))
  
  b = betweenness(g, v=V(g), directed=FALSE)
  show("Betweenness of Non-AP Genes:")
  show(b[!idx])
  show("Betweenness of AP Genes:")
  show(b[idx])
  show("Avg Non-AP Betweenness:")
  show(mean(b[!idx]))
  show("Avg AP Betweenness:")
  show(mean(b[idx]))
}

# plots the network using igraph library
# TODO: look at pathview
# color nodes differently
plot_graph = function(am){
  library(igraph)
  g = graph.adjacency(am)
  plot(g, vertex.size=0, edge.arrow.size=0)
}


