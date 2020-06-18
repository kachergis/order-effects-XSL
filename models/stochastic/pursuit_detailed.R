# Stevens et al. 2014 pursuit model
# The pursuit model assumes that when word w is heard, it’s strongest associated meaning (h) is selected, 
# and the association A(w,h) is updated based on the presence or absence of that referent. 
# A(w,h) is increased if h is present in the current context, and no other association involving w is strengthened.
# If h is not present, A(w,h) is decreased and a single new referent is chosen from the context. 
# Whenever the P(h’|w) for some h’ exceeds a threshold value, that meaning is added to the lexicon (never to be removed).

# for each word w in utterance
# 1) if w is novel, choose a referent in context with minimum association and give it gamma
# 2) else, add new objects to w's hypothesis space
# 3) if w is in the lexicon: if w's hypothesized meaning is in the scene, reward that meaning
#	 otherwise, penalize that meaning
# 4) if w is not in the current lexicon:
#      for each obj o in the scene, if o is not in the current lexicon, reward that meaning
# 5) update the current lexicon (threshold)

# following Bush and Mosteller (1951) Linear Reward-Penalty (LR-P) scheme:
# reward: p(h) = p(h) + gamma*(1-p(h)) and reduce all h'!=h: p(h') = p(h')*(1-gamma)
# penalty: p(h) = p(h)*(1-gamma)

### where did i get the following??
#	for all h'!=h, p(h') = gamma/(n-1) + p(h')*(1-gamma) where n=# of hypotheses being considered

#ord = matrix(c(1,2, 1,3), nrow=2, ncol=2, byrow=T)
#ord = matrix(c(1,2,3, 1,4,5, 2,3,4, 5,6,1), nrow=4, ncol=3, byrow=T)
#params = c(.2, .6, .05)

model <- function(params, ord=c(), ord_name="", reps=1, K=NA, verbose=F) {

	gamma = params[1] # learning rate
	thresh = params[2] # threshold (prob) to move an association to the known lexicon 
	lambda = params[3] # smoothing prob

	voc_sz = max(ord, na.rm=TRUE) # vocabulary size
	ppt = dim(ord)[2] # pairs per trial
	m <- matrix(0, voc_sz, voc_sz) # association matrix A(w,h)
  lexicon <- matrix(0, voc_sz, voc_sz) # for any association that reaches threshold
  compScore = rep(0, reps*nrow(ord))
	freq = rep(0,voc_sz) # number of occurrences per pair, so far (to index the resps matrix)
  
	perf = c()
	for(rep in 1:reps) {
		for(t in 1:dim(ord)[1]) {
			tr = as.integer(ord[t,])
			novel = tr[which(freq[tr]==0)] # novel words
			freq[tr] = freq[tr] + 1 
			have_hypoths = tr[which(rowSums(m[tr,1:voc_sz])!=0)] # words which have non-zero assocs
			#need_hypoths = tr[which(rowSums(m[tr,1:voc_sz])==0)] # words with zero 
			
			# initialize...
			# 1) if w is novel, select one from available refs with minimum strength
			for(w in 1:length(novel)) {
			  max_ref = apply(m[tr,], 1, max) # max association for each referent on trial
			  # select one referent with minimum maximum association
			  min_ref = sample(tr[which(max_ref == min(max_ref))], 1)
				m[novel[w], min_ref] = gamma 
			}
			
			need_hypoths = c() # for any disconfirmed words that need a new referent
			for(w in have_hypoths) {
				hypo = which(m[w,]==max(m[w,1:voc_sz])) # strongest one
				if(length(hypo)>1) hypo = sample(hypo, 1) # (choose one if >1 strongest)
				if(!is.element(hypo, tr)) { # not on trial, so weaken:
					m[w,hypo] = m[w,hypo]*(1-gamma) # disconfirmed
					need_hypoths = c(need_hypoths, w) # not on trial, so choose a new hypothesized referent
				} else {
					m[w,hypo] = m[w,hypo] + gamma*(1-m[w,hypo]) # confirmed: strengthen
					m[w,which(1:voc_sz!=hypo)] = m[w,which(1:voc_sz!=hypo)] * (1-gamma) # weaken others (they do this?)
				}
			}
			store = need_hypoths # reward A(w,h') for a randomly selected new referent (from trial)
			if(length(store)==1) {
			  new_hyps = store
			} else {
			  new_hyps = sample(store, length(store), replace=FALSE)
			}
			for(w in 1:length(store)) {
				m[need_hypoths[w], new_hyps[w]] = gamma*(1-m[need_hypoths[w],new_hyps[w]])
			}
			
			Pm_w = m+lambda
			Pm_w = Pm_w / rowSums(Pm_w) # Eq 1
			lexicon[which(Pm_w>thresh)] = 1 # add to lexicon
			#compScore[t] = sum(diag(lexicon) / rowSums(lexicon+1e-9)) # should use index
			compScore[t] = sum(diag(Pm_w) / rowSums(Pm_w+1e-9))
		}
		#perf = c(perf, sum(diag(lexicon)>0)/voc_sz)
	  perf = c(perf, sum(diag(lexicon) / rowSums(lexicon)))
		if(verbose) print(m)
	}
	#resp_prob = diag(Pm_w) / rowSums(Pm_w)
	resp_prob = diag(lexicon) / rowSums(lexicon)
	want = list(perf=resp_prob, matrix=lexicon, compScore=compScore) # , traj=traj
	return(want)
}


#model(c(.2, .6, .001), orig)

#start_t <- Sys.time()
#simulate_subjects(orig, c(.2, .6, .01))
#end_t <- Sys.time()

#end_t - start_t # 5.2s