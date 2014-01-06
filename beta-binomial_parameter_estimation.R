# method of moments
methodOfMoments <- function(data) {
	size <- nrow(data)
	N <- data[,1]
	K <- data[,2]
	m1 <- sum(K) / size
	m2 <- sum(K^2) / size
	n <- mean(N)
	al <- (n * m1 - m2) / (n * (m2 / m1 - m1 - 1) + m1)
	be <- (n - m1) * (n - m2 / m1) / (n * (m2 / m1 - m1 -1) + m1)
	c(al, be)
}

# maximum likelihood
F1 <- function(data, param) {
	n <- nrow(data) # sample size
	al <- param[1]
	be <- param[2]
	N <- data[,1]
	K <- data[,2]
	f <- n * (digamma(sum(param)) - digamma(al)) + sum(digamma(K + al) - digamma(N + sum(param)))
	g <- n * (digamma(sum(param)) - digamma(be)) + sum(digamma(N - K + be) - digamma(N + sum(param)))
	c(f, g)
}

derivF1 <- function(data, param) { # d/d(param) F
	n <- nrow(data) # sample size
	N <- data[,1]
	K <- data[,2]
	al <- param[1]
	be <- param[2]
	dadf <- n * (trigamma(sum(param)) - trigamma(al)) + sum(trigamma(K + al) - trigamma(N + sum(param)))
	dbdf <- n * trigamma(sum(param)) - sum(trigamma(N + sum(param)))
	dbdg <- n * (trigamma(sum(param)) - trigamma(be)) + sum(trigamma(N - K + be) - trigamma(N + sum(param)))
	matrix(c(dadf, dbdf, dbdf, dbdg), ncol=2, byrow=T)
}

estimate <- function(data, init=c(1,1)) {
	p0 <- init
	for(i in 0:100) {
		p1 <- p0 - solve(derivF1(data, p0)) %*% F1(data, p0)
		if((sum(p1 - p0))^2 < 1e-5) break
		p0 <- p1
		print(i)
	}
	p1
}

# from wiki
males <- seq(0, 12)
families <- c(3, 24, 104, 286, 670, 1033, 1343, 1112, 829, 478, 181, 45, 7)

K <- unlist(sapply(males, function(x) { return(rep(x, families[x+1])) }))
N <- rep(12, length(K))
D <- cbind(N, K)

paramByMM <- methodOfMoments(D)
estimate(D, init=paramByMM)
estimate(D, init=c(1,1))

D <- cbind(rep(40, 5), sapply(rbeta(5, 1, 10), function(x) { rbinom(1, 40, x)}))
D
methodOfMoments(D)
estimate(D, c(1, 1))
estimate(D, methodOfMoments(D))
