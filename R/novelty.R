# Associative Uncertainty- (Novelty) & Familiarity-Biased Model
# George Kachergis  george.kachergis@gmail.com


update_known <- function(m, tr_w, tr_o, startval = .01) {
  tr_assocs = m[tr_w, tr_o]
  tr_assocs[which(tr_assocs==0)] = startval
  m[tr_w, tr_o] = tr_assocs
  
  # for any other experienced word (not on this trial), fill in startval
  
  fam_objs = which(colSums(m)>0)
  fam_words = which(rowSums(m)>0)
  
  for(w in tr_w) {
    zeros = which(m[w,fam_objs]==0)
    m[w,zeros] = startval
  }
  
  for(o in tr_o) {
    zeros = which(m[fam_words,o]==0)
    m[zeros,o] = startval
  }
  
  return(m)
}


model <- function(params, ord=c(), start_matrix=c(), reps=1, test_noise=0) {
	X <- params[1] # associative weight to distribute
	B <- params[2] # weighting of uncertainty vs. familiarity
	C <- params[3] # decay
	
	voc = unique(unlist(ord$words))
	ref = unique(unlist(ord$objs[!is.na(ord$objs)]))
	voc_sz = length(voc) # vocabulary size
	ref_sz = length(ref) # number of objects
	freq_w = rep(0,voc_sz) # freq[i] = times word i has appeared
	freq_o = rep(0,ref_sz)
	names(freq_w) = voc
	names(freq_o) = ref
	traj = list()
	if(is.matrix(start_matrix)) {
	  m <- start_matrix
	} else {
	  m <- matrix(0, voc_sz, ref_sz) # association matrix
	}
	colnames(m) = ref
	rownames(m) = voc
	perf = matrix(0, reps, voc_sz) # a row for each block
	# training
	for(rep in 1:reps) { # for trajectory experiments, train multiple times
	  for(t in 1:length(ord$words)) { 
		#print(format(m, digits=3))
		
	    tr_w = unlist(ord$words[t])
	    tr_w = tr_w[!is.na(tr_w)]
	    tr_w = tr_w[tr_w != ""]
	    tr_o = unlist(ord$objs[t])
	    tr_o = tr_o[!is.na(tr_o)]
		m = update_known(m, tr_w, tr_o) # what's been seen so far?
		
		freq_w[tr_w] = freq_w[tr_w] + 1
		freq_o[tr_o] = freq_o[tr_o] + 1
		novelty_w = 1/(1+freq_w[tr_w])
		novelty_o = 1/(1+freq_o[tr_o])
		
		novelty_w = exp(B*novelty_w)
		novelty_o = exp(B*novelty_o)
		nov = (novelty_w %*% t(novelty_o))
		
		# get all current w,o strengths and normalize to distr X
		assocs = m[tr_w,tr_o]
		denom = sum(assocs * nov)
		m = m*C # decay everything
		
		# update associations on this trial
		m[tr_w,tr_o] = m[tr_w,tr_o] + (X * assocs * (novelty_w %*% t(novelty_o))) / denom 

		index = (rep-1)*length(ord$words) + t # index for learning trajectory
		traj[[index]] = m
		
	  }
	m_test = m+test_noise # test noise constant k
	perf[rep,] = get_perf(m_test)
	}
	
	want = list(perf=perf, matrix=m, traj=traj)
	return(want)
	}

