################################################
# Some usefull functions
stateToId <- function(state){
	# Function stateToId
	# Purpose: Convert state to id
	state <- as.character(state)
	state[state=="B"] <- 1
	state[state=="T"] <- 2
	state[state=="M"] <- 3
	state[state=="R"] <- 4
	state[is.na(state)] <- 0

	return(as.numeric(state))
}

idToState <- function(id){
        # Function idToState
        # Purpose: Convert id to state

    id <- as.numeric(id)
    id[id==1] <- "B"
    id[id==2] <- "T"
    id[id==3] <- "M"
    id[id==4] <- "R"
    id[id==0] <- NA

    return(id)
}
