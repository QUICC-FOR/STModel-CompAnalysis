

M = eq$M
T = eq$T
R = eq$R
B = eq$B
obs = eq$obs

Lik = numeric(length(obs))

Lik[obs=="M"] = M[obs=="M"]
Lik[obs=="B"] = B[obs=="B"]
Lik[obs=="R"] = R[obs=="R"]
Lik[obs=="T"] = T[obs=="T"]

LL = log(Lik)

null = table(obs)/length(obs)
nullLik = numeric(length(obs))
nullLik[obs=="M"] = null[2]
nullLik[obs=="B"] = null[1]
nullLik[obs=="R"] = null[3]
nullLik[obs=="T"] = null[4]
lognullLL = log(nullLik)

1-sum(LL,na.rm=T)/sum(lognullLL,na.rm=T)




