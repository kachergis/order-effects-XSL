# Associative Uncertainty- (Entropy) & Familiarity-Biased Model
# George Kachergis  gkacherg@indiana.edu  June 10, 2011
# modified from filtering model to have remember pairings that appeared on the previous
# trial and separate those from the new ones on this trial (do not cross-associate groups)
# August 5, 2011

shannon.entropy <- function(p) {
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}


update_known <- function(m, tr) {
  startval = .01
  
  for(i in tr) {
    for(c in 1:dim(m)[2]) {
      if(sum(m[,c]>0) & m[i,c]==0) {
        m[i,c] = startval
        m[c,i] = startval
      }
    }
    for(j in tr) {
      if(m[i,j]==0) m[i,j] = startval
      if(m[j,i]==0) m[j,i] = startval
    }
  }
  return(m)
}

model <- function(params, ord=c(), ord_name="", name="model", save_traj=FALSE) {
	X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	alpha <- params[3] # decay
	attnwt <- params[4]
	
	voc_sz = max(unlist(ord), na.rm=T) # vocabulary size
	
	m <- matrix(0, voc_sz, voc_sz) # association matrix
	trial_sz = dim(ord)[2]
	
	# training
	for(t in 1:dim(ord)[1]) { 
		#print(format(m, digits=3))
		tr = unlist(ord[t,])
		
		# temporal contiguity: separate cont & discont
		if(t>1) { 
			contig = intersect(tr, unlist(ord[t-1,])) 
			discontig = setdiff(tr, unlist(ord[t-1,]))
		} else { contig = c(); discontig = tr }
		
		m = update_known(m, tr) # what's been seen so far?
		
		# apportion X (equally, for now--maybe a parm)
		cu = length(contig)^2
		du = length(discontig)^2
		Xc = (X*cu) / (cu+du) 
		Xd = (X*du) / (cu+du)
		
		Xc = Xc*attnwt # to attend more or less to contig pairs
		Xd = Xd*(1-attnwt) # or use alpha as general decay
		
		entC = c() # more entropy = more dispersive
		entD = c()
		for(w in contig) { entC = c(entC, shannon.entropy(m[w,])) }
		for(w in discontig) { entD = c(entD, shannon.entropy(m[w,])) }
		entC = exp(B*entC)
		entD = exp(B*entD)
		
		# get w,o strengths for current trial and normalize to distr X
		discontig <- unlist(discontig)
		contig <- unlist(contig)
		assocsC = m[contig,contig]
		assocsD = m[discontig,discontig]
		
		denomC = sum(assocsC * (entC %*% t(entC)))
		denomD = sum(assocsD * (entD %*% t(entD)))
		
		m = m*alpha # decay everything
		
		m[contig,contig] = m[contig,contig] + (Xc*assocsC*(entC%*%t(entC)))/denomC # update assocs
		m[discontig,discontig] = m[discontig,discontig] + (Xd*assocsD*(entD%*%t(entD)))/denomD
		
		}
	m = m+.01
	#return(test(m, voc_sz, ord_name))
	perf = diag(m) / rowSums(m)
	want = list(perf=perf, matrix=m)
	return(want)
	}

