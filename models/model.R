# Associative Uncertainty- (Entropy) & Familiarity-Biased Model
# George Kachergis  gkacherg@indiana.edu  June 10, 2011

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


model <- function(params, ord, start_matrix=c(), ord_name="", reps=1, name="model", print_matrix=FALSE) {
	X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	C <- params[3] # decay
	
	voc_sz = max(unlist(ord), na.rm=TRUE) # vocabulary size
	traj = list()
	compScore = rep(0, nrow(ord))
	ppt = ncol(ord) # pairs per trial
	mean_ent = c()
	
	
	if(length(start_matrix)==0) { # no prior experience (e.g., pre-trained associations)
	  m <- matrix(0, voc_sz, voc_sz) # empty association matrix
	} else {
	  m = start_matrix # pretrained association
	}
	
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord)) { 
		#print(format(m, digits=3))
		
		tr = as.integer(ord[t,])
		tr = tr[!is.na(tr)]
		m = update_known(m, tr) # what's been seen so far?
		ent = c() # more entropy = more dispersive
		for(w in tr) { ent = c(ent, shannon.entropy(m[w,])) }
		mean_ent = c(mean_ent, sum(ent)/ppt)
		ent = exp(B*ent)
		
		# get all current w,o strengths and normalize to distr X
		assocs = m[tr,tr]
		
		denom = sum(assocs * (ent %*% t(ent)))
		
		m = m*C # decay everything
		
		m[tr,tr] = m[tr,tr] + (X * assocs * (ent %*% t(ent))) / denom # update assocs
		
		traj[[t]] = m
		compScore[t] = sum(diag(m) / rowSums(m+1e-9)) #sum(diag(m))
	  }
	}
	perf = diag(m) / rowSums(m)
	want = list(perf=perf, matrix=m, traj=traj, compScore=compScore)
	return(want)
	}

