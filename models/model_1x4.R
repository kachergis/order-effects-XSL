# Associative Uncertainty- (Entropy) & Familiarity-Biased Model
# George Kachergis  gkacherg@indiana.edu  June 10, 2011

shannon.entropy <- function(p) {
	if (min(p) < 0 || sum(p) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
	}

update_known <- function(m, word, objs) {
	startval = .01
	for(c in 1:dim(m)[2]) {
		if(sum(m[,c]>0) & m[word,c]==0) {
			m[word,c] = startval
		}
	}
	for(j in objs) {
		if(m[word,j]==0) m[word,j] = startval
	}
	return(m)
	}

# special version for 1 word x N object trials (Koehne et al 2013)
# each row: w o o o o
model <- function(params, ord, ord_name="", reps=1, name="model", print_matrix=FALSE) {
	X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	C <- params[3] # decay
	
	voc_sz = max(unlist(ord[,1]), na.rm=TRUE) # vocabulary size
	nrefs = max(unlist(ord), na.rm=TRUE) 
	traj = list()
	compScore = rep(0, nrow(ord))
	ppt = ncol(ord)-1 # pairs per trial
	mean_ent = c()
	m <- matrix(0, voc_sz, nrefs) # association matrix
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord)) { 
		#print(format(m, digits=3))
		
	  word = as.integer(ord[t,1])
		objs = as.integer(ord[t,2:ncol(ord)])
		m = update_known(m, word, objs) # what's been seen so far?
		ent = c() # more entropy = more dispersive
		for(o in objs) { ent = c(ent, shannon.entropy(m[,o])) }
		mean_ent = c(mean_ent, sum(ent)/ppt)
		ent = exp(B*ent)
		
		# get all current w,o strengths and normalize to distr X
		assocs = m[word,objs]
		
		denom = sum(assocs * ent)
		
		m = m*C # decay everything
		
		m[word,objs] = m[word,objs] + (X * assocs * ent) / denom # update assocs
		
		traj[[t]] = m
		compScore[t] = sum(diag(m) / rowSums(m+1e-9)) #sum(diag(m))
	  }
	}
	perf = diag(m) / rowSums(m)
	want = list(perf=perf, matrix=m, traj=traj, compScore=compScore)
	return(want)
	}

