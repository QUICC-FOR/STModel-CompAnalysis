# optimisation
# prerun the model with glm to find some parameter values that are not jointly constrained
# use logits

model = function(params, dat, step= 1) 

{
    st0 = dat$st0
    st1 = dat$st1
    ET = dat$ET
    EB = dat$EB
    EM = dat$EM
    ENV1 = dat$ENV1
    ENV2 = dat$ENV2
    itime = dat$itime
    
    names(params) = c("ab0", "ab1", "ab2", "ab3","ab4",
"at0", "at1" , "at2", "at3", "at4",
"bb0", "bb1", "bb2", "bb3", "bb4", 
"bt0", "bt1", "bt2", "bt3", "bt4", 
"tt0", "tt1", "tt2", "tt3", "tt4", 
"th0", "th1", "th2", "th3", "th4", 
"e0", "e1", "e2", "e3", "e4")

    
	lik = numeric(length(st0))

    alphab 	= params["ab0"] + params["ab1"]*ENV1 + params["ab2"]*ENV2 + params["ab3"]*ENV1^2 + params["ab4"]*ENV2^2 
    alphat 	= params["at0"] + params["at1"]*ENV1 + params["at2"]*ENV2 + params["at3"]*ENV1^2 + params["at4"]*ENV2^2 
    betab 	= params["bb0"] + params["bb1"]*ENV1 + params["bb2"]*ENV2 + params["bb3"]*ENV1^2 + params["bb4"]*ENV2^2 
    betat 	= params["bt0"] + params["bt1"]*ENV1 + params["bt2"]*ENV2 + params["bt3"]*ENV1^2 + params["bt4"]*ENV2^2 
    theta	= params["th0"] + params["th1"]*ENV1 + params["th2"]*ENV2 + params["th3"]*ENV1^2 + params["th4"]*ENV2^2 
    thetat	= params["tt0"] + params["tt1"]*ENV1 + params["tt2"]*ENV2 + params["tt3"]*ENV1^2 + params["tt4"]*ENV2^2 
    eps 	= params["e0"]  + params["e1"]*ENV1 + params["e2"]*ENV2  + params["e3"]*ENV1^2 + params["e4"]*ENV2^2 
    
	# Compute the likelihood of observations
	pBM = (betat*(ET+EM)*(1-eps))
	pBR = eps
	pBB = (1 - eps - betat*(ET+EM)*(1-eps))

	pTT = (1 - eps - betab*(EB+EM)*(1-eps))
	pTM = (betab*(EB+EM)*(1-eps))
	pTR = eps
	
	pMB = (theta*(1-thetat)*(1-eps))
	pMT= (theta*thetat*(1-eps))

	pMM = ((1 - eps)*(1 - theta))

	pMR = eps
	phib = alphab*(EM + EB)*(1-alphat*(ET+EM))
	phit = alphat*(EM + ET)*(1-alphab*(EB+EM))
	phim = alphab*(EM + EB)*alphat*(EM + ET)
	lik[st0 == "R" & st1 == "B"] = phib[st0 == "R" & st1 == "B"] 	
	lik[st0 == "R" & st1 == "T"] = phit[st0 == "R" & st1 == "T"]	
	lik[st0 == "R" & st1 == "M"] = phim[st0 == "R" & st1 == "M"] 			
	lik[st0 == "R" & st1 == "R"] = (1 - phib - phit - phim)[st0 == "R" & st1 == "R"] 


    # lik might be <0 !!!!
    # for T->T, M->M, B->B and R->R transitions
    # it will give NaN when log
    

    # lik might be equal = 0 -> give -Inf at log   
    # for instance when neighbor (seeds) = 0
    #lik[lik == 0] = .Machine$double.xmin
    
    # calculate sum log likelihood
    # lik might be <0 (in the 1 - other proba)
#    if(sum(lik<0)>0)
#    {
#    sumLL = -.Machine$double.xmax
#    }else{
    sumLL = sum(log(lik))
#    }
    
    if(is.infinite(sumLL)) sumLL = -.Machine$double.xmax
	if(is.nan(sumLL)) sumLL = -.Machine$double.xmax
	if(is.na(sumLL)) print("sumLL: na values!")
	
#	print(sumLL)
    # return a value to minimize in GenSA
	return(-sumLL)
}


#-----------------------------------------------------------------------------

