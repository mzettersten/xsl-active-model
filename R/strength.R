# Associative Familiarity-Biased Model
# (just the familiarity/strength bias of the Kachergis et al. 2012 model)
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


model <- function(params, ord=c(), reps=1, test_noise=0) {
	X <- params[1] # associative weight to distribute
	C <- params[2] # decay
	
	voc = unique(unlist(ord$words))
	ref = unique(unlist(ord$objs[!is.na(ord$objs)]))
	voc_sz = length(voc) # vocabulary size
	ref_sz = length(ref) # number of objects
	traj = list()
	m <- matrix(0, voc_sz, ref_sz) # association matrix
	perf = matrix(0, reps, voc_sz) # a row for each block
	
	colnames(m) = ref
	rownames(m) = voc
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
		
		# get all current w,o strengths and normalize to distr X
		assocs = m[tr_w,tr_o]
		denom = sum(assocs)
		m = m*C # decay everything
		# update associations on this trial
		m[tr_w,tr_o] = m[tr_w,tr_o] + (X * assocs) / denom 

		index = (rep-1)*length(ord$words) + t # index for learning trajectory
		traj[[index]] = m
	  }
	m_test = m+test_noise # test noise constant k
	perf[rep,] = get_perf(m_test)
	}
	want = list(perf=perf, matrix=m, traj=traj)
	return(want)
	}

