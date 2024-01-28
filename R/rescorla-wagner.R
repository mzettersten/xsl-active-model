# Rescorla-Wagner (1972) associative learning model
# adapted for cross-situational word learning by
# George Kachergis  george.kachergis@gmail.com

model <- function(params, ord=c(), start_matrix=c(), reps=1, test_noise=0) {
	C = params[2] # decay
	lambda = params[3] # maximum associative value that a CS can achieve - should be larger than learning rate
	beta = params[1]*lambda # learning rate -- a proportion of lambda
	alpha = 1 # salience (fix at 1 unless manipulated)

	voc_sz = max(unlist(ord$words), na.rm=TRUE) # vocabulary size
	ref_sz = max(unlist(ord$objs), na.rm=TRUE) # number of objects
	traj = list()
	if(is.matrix(start_matrix)) {
	  m <- start_matrix
	} else {
	  m <- matrix(0, voc_sz, ref_sz) # association matrix
	}
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

		# if objects are cues that predict words, then we want colSums;
		# if words are cues, use rowSums--but only of the currently-presented stimuli
		if(length(tr_w)==1) {
		  pred = m[tr_w,tr_o] # should prediction be based on entire col of obj assocs?
		} else if(length(tr_o)==1) {
		  pred = sum(m[tr_w,tr_o])
		} else {
		  pred = colSums(m[tr_w,tr_o])
		}
		delta = alpha*beta*(lambda - pred)
		m[tr_w,tr_o] = m[tr_w,tr_o] + delta

		m = m*C

		index = (rep-1)*nrow(ord$words) + t # index for learning trajectory
		traj[[index]] = m
	  }

	m_test = m+test_noise # test noise constant k
	perf[rep,] = diag(m_test) / rowSums(m_test)
	}
	want = list(perf=perf, matrix=m, traj=traj)
	return(want)
	}

