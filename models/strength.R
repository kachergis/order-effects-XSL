# Associative Familiarity-Biased Model
# (just the familiarity/strength bias of the Kachergis et al. 2012 model)
# George Kachergis  george.kachergis@gmail.com

update_known <- function(m, tr_w, tr_o) {
  startval = .01
  
  for(i in tr_w) {
    for(c in 1:dim(m)[2]) {
      if(sum(m[,c]>0) & m[i,c]==0) {
        m[i,c] = startval
        if(c<nrow(m)) m[c,i] = startval
      }
    }
    for(j in tr_o) {
      if(m[i,j]==0) m[i,j] = startval
      if(j<nrow(m) && m[j,i]==0) m[j,i] = startval
    }
  }
  return(m)
}


model <- function(params, ord=c(), reps=1, test_noise=.01) {
	X <- params[1] # associative weight to distribute
	C <- params[2] # decay
	
	voc_sz = max(unlist(ord), na.rm=TRUE) # vocabulary size
	ref_sz = voc_sz # number of objects
	ppt = ncol(ord) # pairs per trial
	traj = list()
	compScore = rep(0, nrow(ord)*reps)
	m <- matrix(0, voc_sz, ref_sz) # association matrix
	perf = matrix(0, reps, voc_sz) # a row for each block
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:nrow(ord)) { 
		#print(format(m, digits=3))
		
		tr_w = as.integer(ord[t,])
		tr_w = tr_w[!is.na(tr_w)]
		tr_o = as.integer(ord[t,])
		tr_o = tr_o[!is.na(tr_o)]
		m = update_known(m, tr_w, tr_o) # what's been seen so far?
		
		# get all current w,o strengths and normalize to distr X
		assocs = m[tr_w,tr_o]
		denom = sum(assocs)
		m = m*C # decay everything
		# update associations on this trial
		m[tr_w,tr_o] = m[tr_w,tr_o] + (X * assocs) / denom 

		index = (rep-1)*length(ord$trials) + t # index for learning trajectory
		traj[[index]] = m
		compScore[index] = sum(diag(m) / rowSums(m+1e-9))
	  }
	m_test = m+test_noise # test noise constant k
	perf[rep,] = diag(m_test) / rowSums(m_test)
	}
	want = list(perf=perf[reps,], matrix=m, traj=traj, compScore=compScore)
	return(want)
	}

