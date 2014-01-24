# code to do network analysis of AP 
# genes in the context of Reactome.bd
# networks.
# Under construction.

# load libraries
library(hgu133plus2.db)
library(antiProfiles)
library(plyr)
library(glasso)
library(reactome.db)
library(RCytoscape)
library(ppcor)

# loads data into the workspace as global vars 
# c_exprs and c_condition
get_colon_data = function(){
  #library(hgu133plus2.db)
  load("apData_gpl570.rda")
  idx = pData(apData_gpl570)$Tissue == "colon"
  c_exprs <<- exprs(apData_gpl570)[,idx]
  #c_condition <<- pData(apData_gpl570)[idx,]$Status
  c_condition <<- as.numeric(pData(apData_gpl570)[idx,]$SubType == "colorectal_cancer")
}

# Takes expression data in, indexed by probe_ids
# then converts them to Entrez ids, 
# collapsing the multiple probes that match to one
# entrez ids via the median
convert_probes_entrez = function(expr_data){
  #library(hgu133plus2.db)
  
  probe_ids = rownames(expr_data)
  entrez_ids = hgu133plus2ENTREZID[probe_ids]
  entrez_ids = toTable(entrez_ids)
  
  #library(plyr)
  
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
  #library(hgu133plus2.db)
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
  #library(hgu133plus2.db)
  
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

# prints gets several pathways, and runs glasso on them
# prints the results to a chart
sample_paths_glasso = function(){
  #library(antiProfiles)
  #library(glasso)
  #library(hgu133plus2.db)
  #library(reactome.db)
  
  #get_colon_data()	# should use setup_colon_ap_data instead
  
  # pathways list:
  # apoptosis, degredation of the extracellular matrix, growth hormone receptor signaling (immune -> cytokine),
  # signaling by EGFR, Signalling by NOTCH1
  sample_paths = c("169911", "1474228", "982772", "177929", "2644602")
  
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


fit_glasso_to_path = function(pathid){
  control_samples = c_condition == 0
  cancer_samples = c_condition == 1
  #data = convert_probes_entrez(data)
  data = load_path_exprs_subset(pathid)
  data = convert_probes_entrez(data)
  r = rel_from_pathid_slim2(pathid, data)
  rel = r[[1]]
  
  zero_partial_corrs = r[[2]]
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
    
  controlg = glasso(cov_control, rho=0.001, zero=zero_partial_corrs)    
  cancerg = glasso(cov_cancer, rho=0.001, zero=zero_partial_corrs)
  show("glasso done")

  #return(c(controlg, cancerg, dapgl))
  return(list(controlg, cancerg, dapgl))
}

fit_pcor_to_path = function(pathid){
  control_samples = c_condition == 0
  cancer_samples = c_condition == 1

  data = load_path_exprs_subset(pathid)
  data = convert_probes_entrez(data)

  dapgl = ap_entrez_ids %in% rownames(data)
  dapgl = ap_entrez_ids[dapgl]

  
  am = get_graph_am(pathid)
  #con_pc = get_partial_correlations(data[,control_samples], am)
  can_pc = get_partial_correlations(data[,cancer_samples], am)
  
  return(list(con_pc, can_pc, dapgl, am, data))
}

fit_pcor_lm_to_path = function(pathid){
  control_samples = c_condition == 0
  cancer_samples = c_condition == 1

  data = load_path_exprs_subset(pathid)
  data = convert_probes_entrez(data)

  dapgl = ap_entrez_ids %in% rownames(data)
  dapgl = ap_entrez_ids[dapgl]

  
  am = get_graph_am(pathid)
  con_pc = get_partial_correlations_lm(data[,control_samples], am)
  can_pc = get_partial_correlations_lm(data[,cancer_samples], am)
  
  return(list(con_pc, can_pc, dapgl, am, data))
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
  sample_paths <<- c("169911", "1474228", "982772", "177929", "2644602")
  # sample paths found by manually querying the file "reactome_path_ids.txt"
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
# note: igraph doesn't seem to plot nicely (at least not on my computer)
plot_graph = function(am){
  library(igraph)
  g = graph.adjacency(am)
  plot(g, vertex.size=1, edge.arrow.size=.1, layout=layout.fruchterman.reingold)
  #plot(g)
}


# has problems with taking the inverse of a singular matrix to find the partial
# covariance when using the pcor.test method
get_partial_correlations = function(exprs, am){
  pcm = am #partial correlation matrix
  pcm[lower.tri(pcm)] = 0 # make lower triangle 0 to prevent work duplication
  z = 0
  for (i in 1:length(rownames(am))){
    print("==================================")
    neighbors = which(pcm[i,]!=0, arr.ind=T)
    eye = diag(length(neighbors))
    #print(neighbors)
    if (length(neighbors) == 0){
      #print("skip")
      next
    }

    if (length(neighbors) > 2){
      #print("many")
      #print(neighbors)
      for (j in 1:length(neighbors)){
        #print(i)
        #print(neighbors[j])
        r = try(pcor.test(t(exprs[i,]), t(exprs[neighbors[j],]), t(exprs[neighbors[1-eye[j,] == 1],]), method="spearman"))
        #print(r)
        if (class(r) == "try-error"){
          print("error_caught")
          z = z + 1
          pcm[i,neighbors[j]] = 0
        } else {
          pcm[i,neighbors[j]] = r$estimate
        }
      }

    } else if (length(neighbors) == 2) {
      print("two")
      r = pcor(t(exprs[append(i, neighbors),]))
      pcm[i, neighbors[1]] = r$estimate[1,2]
      pcm[i, neighbors[2]] = r$estimate[1,3]
    
    } else {
      print("one")
      #print(i)
      #print(neighbors[1])
      r = cor(exprs[i,], exprs[neighbors[1],])
      pcm[i, neighbors[1]] = r
    }
  }
  print(z)
  return(pcm)
}


get_partial_correlations_lm = function(exprs, am){
  pcm = am #partial correlation matrix
  pcm[lower.tri(pcm)] = 0 # make lower triangle 0 to prevent work duplication

  dfexprs = as.data.frame(t(exprs))
  
  z = 0
  for (i in 1:length(rownames(am))){
    print("==================================")
    
    name_i = rownames(am)[i]
    links = which(pcm[i,]!=0, arr.ind=T)
    
    if (length(links) == 0){
      next
    }
    
    for (j in 1:length(links)){
      name_j = rownames(am)[links[j]]
      
      print("---")
      print(name_i)
      print(name_j)
      print("===")
      
      i_neighbors = am[i,]
      i_neighbors[j] = 0
      i_neighbors = rownames(am)[i_neighbors==1]
      
      j_neighbors = am[j,]
      j_neighbors[i] = 0
      j_neighbors = rownames(am)[j_neighbors==1]
      
      i_residuals = dfexprs[name_i]
      j_residuals = dfexprs[name_j]
      
      if (length(i_neighbors) > 0){
        i_form = paste('`', name_i, '` ~ `', paste(i_neighbors, sep='', collapse='` + `'), '`', sep='')
        i_m = lm(i_form, dfexprs)
        i_residuals = i_m$residuals
      }
      
      if (length(j_neighbors) > 0){
        j_form = paste('`', name_j, '` ~ `', paste(j_neighbors, sep='', collapse='` + `'), '`', sep='')
        j_m = lm(j_form, dfexprs)
        j_residuals = j_m$residuals
      }
      
      pcm[name_i, name_j] = cor(i_residuals, j_residuals)[1]
      
    }
  }  
  return(pcm)
}


plot_cytoscape = function(am, controlg, cancerg, ap_ids){
  #library(RCytoscape)
  g <- as(am, "graphNEL")
  g <- initEdgeAttribute(g, "weight", "numeric", 1.0)
  #try(cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=TRUE))
  cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=TRUE)
  #g <- initEdgeAttribute(g, "weight", "numeric", 1.0)
  #cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=FALSE)

  print("Base Graph Sent")
  displayGraph(cw)

  layoutNetwork(cw, layout.name="force-directed")
  # check getLayoutNames(cw) for other options
  redraw(cw)

  print("Graph layed-out")

  # how to set the color of the nodes...
  #> setNodeAttributesDirect(cw, "gft", "numeric", rownames(a), 1*(rownames(a)> "5000"))
  #> gftcolors = c('#00AA00', '#00FF00')
  #> setNodeColorRule(cw, 'gft', c(0, 1), gftcolors, mode='interpolate')
  #redraw(cw)


  # Ideal graph:
  # Nodes are squares or circles if the are ap or not
  # nodes colored by delta of expression levels
  # links are colored based on the delta of strength of the positive or negative correlation
  # second graph
  # simpler
  # still squares and circles for nodes via ap or not
  # node color based on expression level
  # link color based on strength of positive or negative correlation
  
  #i = which(am!=0, arr.ind=T)
  #am[i] # gets values
  
  # get indicies
  a = am
  a[lower.tri(a)] = 0 # gets rid of idicies where row>col
  i = which(a!=0, arr.ind=T) # gets indicies
  i = i[order(i[,1], i[,2]),] # sorts by row ascending order
  
  edge_names = getAllEdges(cw)
  edge_names = rev(edge_names) # sorted by row ascending

  # set the background color to white
  setDefaultBackgroundColor(cw, "#FFFFFF", "default")
  
  # set node shapes
  #print(ap_ids)
  #print(getNodeShapes(cw))
  setNodeShapeDirect(cw, ap_ids, rep('diamond', length(ap_ids)))
  
  print("AP genes made diamond")
  
  # Set edge colors
  #setEdgeColorDirect(cw, edge_names, "#FF0000")
  #setEdgeColorDirect(cw, edge_names[controlg$wi[i] > 0], "#00FF00")
  #setEdgeColorDirect(cw, edge_names[controlg$wi[i] == 0], "#000000")
  #redraw(cw)
  
  setEdgeAttributesDirect(cw, "strength", "numeric", edge_names, round(10*controlg$wi[i]))
  control_points = c(-20.0, 0.0, 20.0); # TODO: These are probably not the best points
  colors = c("#FF0000", "#000000", "#00FF00");
  setEdgeColorRule(cw, "strength", control_points, colors, mode="interpolate")
  #setEdgeColorRule(cw, "strength", c(-1, 1), c("#FF0000", "#00FF00"), mode="interpolate")
  
  displayGraph(cw)
  redraw(cw)
  print("Edges Colored")
  
  
  # set edge widths
  #setEdgeLineWidthDirect(cw, edge_names, abs(controlg$wi[i]))
      
  print("Edge Widths Set")
  
}


plot_cy_single = function(am, glasso_out, ap_ids){
  #library(RCytoscape)
  g <- as(am, "graphNEL")
  g <- initEdgeAttribute(g, "weight", "numeric", 1.0)
  #try(cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=TRUE))
  cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=TRUE)
  #g <- initEdgeAttribute(g, "weight", "numeric", 1.0)
  #cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=FALSE)

  print("Base Graph Sent")
  displayGraph(cw)

  layoutNetwork(cw, layout.name="force-directed")
  # check getLayoutNames(cw) for other options
  redraw(cw)

  print("Graph layed-out")

  # how to set the color of the nodes...
  #> setNodeAttributesDirect(cw, "gft", "numeric", rownames(a), 1*(rownames(a)> "5000"))
  #> gftcolors = c('#00AA00', '#00FF00')
  #> setNodeColorRule(cw, 'gft', c(0, 1), gftcolors, mode='interpolate')
  #redraw(cw)


  # Ideal graph:
  # Nodes are squares or circles if the are ap or not
  # nodes colored by delta of expression levels
  # links are colored based on the delta of strength of the positive or negative correlation
  # second graph
  # simpler
  # still squares and circles for nodes via ap or not
  # node color based on expression level
  # link color based on strength of positive or negative correlation
  
  #i = which(am!=0, arr.ind=T)
  #am[i] # gets values
  
  # get indicies
  a = am
  a[lower.tri(a)] = 0 # gets rid of idicies where row>col
  i = which(a!=0, arr.ind=T) # gets indicies
  i = i[order(i[,1], i[,2]),] # sorts by row ascending order
  
  edge_names = getAllEdges(cw)
  edge_names = rev(edge_names) # sorted by row ascending

  # set the background color to white
  setDefaultBackgroundColor(cw, "#FFFFFF", "default")
  
  # set node shapes
  #print(ap_ids)
  #print(getNodeShapes(cw))
  setNodeShapeDirect(cw, ap_ids, rep('diamond', length(ap_ids)))
  
  print("AP genes made diamond")
  
  # Set edge colors
  #setEdgeColorDirect(cw, edge_names, "#FF0000")
  #setEdgeColorDirect(cw, edge_names[controlg$wi[i] > 0], "#00FF00")
  #setEdgeColorDirect(cw, edge_names[controlg$wi[i] == 0], "#000000")
  #redraw(cw)
  
  setEdgeAttributesDirect(cw, "strength", "numeric", edge_names, round(10*glasso_out$wi[i]))
  control_points = c(-20.0, 0.0, 20.0); # TODO: These are probably not the best points
  colors = c("#FF0000", "#000000", "#00FF00");
  setEdgeColorRule(cw, "strength", control_points, colors, mode="interpolate")
  #setEdgeColorRule(cw, "strength", c(-1, 1), c("#FF0000", "#00FF00"), mode="interpolate")
  
  displayGraph(cw)
  redraw(cw)
  print("Edges Colored")
  
  
  # set edge widths
  #setEdgeLineWidthDirect(cw, edge_names, abs(controlg$wi[i]))
      
  print("Edge Widths Set")
  
}

plot_cy_pcor = function(am, pcm, ap_ids){
  #library(RCytoscape)
  g <- as(am, "graphNEL")
  g <- initEdgeAttribute(g, "weight", "numeric", 1.0)
  #try(cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=TRUE))
  cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=TRUE)
  #g <- initEdgeAttribute(g, "weight", "numeric", 1.0)
  #cw <- new.CytoscapeWindow('ap_net', graph=g, overwriteWindow=FALSE)

  print("Base Graph Sent")
  displayGraph(cw)

  layoutNetwork(cw, layout.name="force-directed")
  # check getLayoutNames(cw) for other options
  redraw(cw)

  print("Graph layed-out")

  # how to set the color of the nodes...
  #> setNodeAttributesDirect(cw, "gft", "numeric", rownames(a), 1*(rownames(a)> "5000"))
  #> gftcolors = c('#00AA00', '#00FF00')
  #> setNodeColorRule(cw, 'gft', c(0, 1), gftcolors, mode='interpolate')
  #redraw(cw)


  # Ideal graph:
  # Nodes are squares or circles if the are ap or not
  # nodes colored by delta of expression levels
  # links are colored based on the delta of strength of the positive or negative correlation
  # second graph
  # simpler
  # still squares and circles for nodes via ap or not
  # node color based on expression level
  # link color based on strength of positive or negative correlation
  
  #i = which(am!=0, arr.ind=T)
  #am[i] # gets values
  
  # get indicies
  a = am
  a[lower.tri(a)] = 0 # gets rid of idicies where row>col
  i = which(a!=0, arr.ind=T) # gets indicies
  i = i[order(i[,1], i[,2]),] # sorts by row ascending order
  
  edge_names = getAllEdges(cw)
  edge_names = rev(edge_names) # sorted by row ascending

  # set the background color to white
  setDefaultBackgroundColor(cw, "#FFFFFF", "default")
  
  # set node shapes
  #print(ap_ids)
  #print(getNodeShapes(cw))
  setNodeShapeDirect(cw, ap_ids, rep('diamond', length(ap_ids)))
  
  print("AP genes made diamond")
  
  # Set edge colors
  #setEdgeColorDirect(cw, edge_names, "#FF0000")
  #setEdgeColorDirect(cw, edge_names[controlg$wi[i] > 0], "#00FF00")
  #setEdgeColorDirect(cw, edge_names[controlg$wi[i] == 0], "#000000")
  #redraw(cw)
  
  setEdgeAttributesDirect(cw, "strength", "numeric", edge_names, pcm[i])
  control_points = c(-.10, 0.0, .10); # TODO: These are probably not the best points
  colors = c("#FF0000", "#000000", "#00FF00");
  setEdgeColorRule(cw, "strength", control_points, colors, mode="interpolate")
  #setEdgeColorRule(cw, "strength", c(-1, 1), c("#FF0000", "#00FF00"), mode="interpolate")
  
  displayGraph(cw)
  redraw(cw)
  print("Edges Colored")
  
  
  # set edge widths
  #setEdgeLineWidthDirect(cw, edge_names, abs(controlg$wi[i]))
      
  print("Edge Widths Set")
  
}


 
