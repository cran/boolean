require("methods")
require("stats")

setClass("booltest",representation(Calculus="character", LogLik="numeric",Variables="vector",Coefficients="vector", StandardErrors="vector", Iterations="numeric", Hessian = "matrix", Gradient = "vector", Zscore = "vector", Probz = "vector", Conf95lo = "vector", Conf95hi = "vector", pstructure = "character",method="character")) 

setMethod("summary", "booltest", function (object) { 
      	z0 <- paste("Model: ", object@pstructure) 
      	z1 <- data.frame(LogLik = object@LogLik, Iterations = object@Iterations) 
      	z2 <- data.frame( Vars=object@Variables,Coefs = round(object@Coefficients, 5), StdErrs = round(object@StandardErrors, 5), Zscore = round(object@Zscore, 3), 
          ProbZ = round(object@Probz, digits = 4),CI95Lo = round(object@Conf95lo, 5), CI95Hi = round(object@Conf95hi, 5)) 
      	format(z2) 
      	cat(z0, "\n") 
      	cat("   \n") 
      	print(z1) 
      	cat("   \n") 
      	print(z2) 
} 
) 


setMethod("plot", "booltest", function(x, y=0, panel="boolfirst") {
         
    	if((panel!="boolfirst")&(panel!="boolprof")){ 
      	stop("Please enter boolfirst or boolprof as panel type.")} 
    	varvec <- x@Variables 
    	if(panel=="boolfirst"){varvec<-subset(varvec, varvec!="cons")} 
    	inst <- rep(0, length(varvec))  # Produce vector of instance numbers 
    	t<-2 
    	v<-1 
    	while(t<length(varvec)+1) { 
      		while(v<t){ 
        		if((varvec[t]==varvec[v])&(inst[v]!=0)){inst[t]<-inst[v]+1} 
        		if((varvec[t]==varvec[v])&(inst[v]==0)){inst[v]<-1;inst[t]<-2} 
        		v<-v+1
		} 
      		t<-t+1 
      		v<-1
	} 
    	gphdim <- ceiling(sqrt(length(varvec)))  # Set graph dimensions 
    	par(mfrow=c(gphdim,gphdim), pty="s") 
    	t<-1                                       # Graph all panels 
    	while(t<length(varvec)+1) { 
        	gphcmd<-paste(panel, "(x, varvec[t], instance=inst[t])", sep="") 
      		eval(parse(text=gphcmd)) 
      		t<-t+1
	} 
} 
) 

setMethod("show", "booltest", function (object) { 
      
	z0 <- paste("Model: ", object@pstructure) 
	z1 <- data.frame(LogLik = object@LogLik, Iterations = object@Iterations) 
      	z2 <- data.frame( Vars=object@Variables, Coefs = round(object@Coefficients, 5), StdErrs = round(object@StandardErrors, 5), 
          Zscore = round(object@Zscore, 3), ProbZ = round(object@Probz, digits = 4), CI95Lo = round(object@Conf95lo, 5), 
          CI95Hi = round(object@Conf95hi, 5)) 
      	z3 <- rbind(object@Variables, round(object@Gradient, 9)) 
      	format(z2) 
      	cat(z0, "\n") 
      	cat("   \n") 
      	print(z1) 
      	cat("   \n") 
      	print(z2) 
      	cat("   \n") 
      	cat("Hessian:\n") 
      	print(object@Hessian) 
      	cat("   \n") 
      	cat("Gradient:\n") 
      	print(z3) 
} 
) 

setMethod("coef", "booltest", function (object) { 
      
	cvec <- data.frame(rbind(object@Variables, 
	round(object@Coefficients, 5))) 
     	 format(cvec) 
} 
) 


booltest <- function(calc, loglik, vars, coefs, ses, iter, hess, grad, zsc, pz, c95lo, c95hi, pstruc,meth)  {

	x<-new("booltest", Calculus=calc, LogLik=loglik, Variables=vars, Coefficients=coefs, StandardErrors=ses, Iterations=iter, Hessian = hess, Gradient = grad, Zscore = zsc, Probz = pz, Conf95lo = c95lo, Conf95hi = c95hi, pstructure = pstruc, method = meth) 
   	return(x) 
} 



boolean <- function(structure,method,maxoptions = "",optimizer="nlm",safety=1,bootstrap=FALSE,bootsize=100,popsize=5000) {

	#checking to ensure that input is of the right form

	if(!is.character(structure)) {structure <- deparse(structure)}
	structure <- gsub(" ","",structure,extended=FALSE)
	if(!is.character(method)) {cat(paste("`method' must be of type 'character' \n"));break}
	if(!is.character(maxoptions)) {cat(paste("`maxoptions' must be of type 'character' \n"));break}
	if ((method != "logit")&(method != "probit")) {cat(paste("Error:  `method' must be Logit or Probit \n"));break}

	#translating input of form y ~ (...) into the dependent variable and probability structure, kept as strings `devpar' and `structure'
	pstructure <- structure
	t <- as.integer(paste(regexpr("~",structure,extended=FALSE)))
	temp <- substr(structure,1,(t-1))
	depvar <- sub(" ","",temp,extended=FALSE)
	structure <- as.character(substring(structure,(t+1),nchar(structure)))

	# correcting for quirk in deparse - string of length n is compressed to a single string

	tempindex <- length(structure)
	holder <- "qghgh"
	for (i in 1:length(structure)) {holder <- sub("qghgh",paste(structure[i],"qghgh",sep=""),holder,extended=FALSE)}
	structure <- sub("qghgh","",holder,extended=FALSE)
	holder <- "qghgh"
	for (i in 1:length(depvar)) {holder <- sub("qghgh",paste(depvar[i],"qghgh",sep=""),holder,extended=FALSE)}
	depvar <- sub("qghgh","",holder,extended=FALSE)
	
	# from probability structure, the variable names are deduced and stored to `vars'

	temp <- gsub(")","",structure,extended=FALSE)
	temp <- gsub("(","",temp,extended=FALSE)
	temp <- gsub("+","&",temp,extended=FALSE)
	temp <- gsub("|","&",temp,extended=FALSE)
	temp <- gsub(" ","",temp,extended=FALSE)
	vars <- character()

	start <- 1
	i <- 1
	t <- as.integer(paste(regexpr("&",temp,extended=FALSE)))

	while (start < nchar(temp)) { 
	
		if (t < 0) {t <- nchar(temp)}
		vars[i] <- substr(temp,start,t)
		temp <- sub("&","x",temp,extended=FALSE)
		i <- i+1
		start <- t+1
		t <- as.integer(paste(regexpr("&",temp,extended=FALSE)))
	
	}


	vars <- sub("&","",vars,extended=FALSE)


	# data set reduced (temporarily) to exclude vectors with missing values

	natemp <- vector()
	eval(parse(text=paste("nahelp <- !is.na(",depvar,")",sep=""))) 
	for (i in 1:length(vars)) {
		if (vars[i] != "cons") {eval(parse(text=paste("natemp <- ",vars[i],sep="")))
		nahelp <- nahelp*!is.na(natemp)}
	
	}

	nahelp <- which(as.logical(nahelp))
	modelvarsnum <- 0

	for (i in 1:length(vars)) {
		if (vars[i] != "cons") {
			if (sum(vars[i] == vars[1:i])==1){
				modelvarsnum <- modelvarsnum + 1
				j <- modelvarsnum
				eval(parse(text=paste("natemp <- ",vars[i],sep="")))
				eval(parse(text=paste("holds",j," <- vector()",sep="")))
				eval(parse(text=paste("holds",j," <- natemp",sep="")))
				natemp <- natemp[nahelp]
				eval(parse(text=paste(vars[i]," <- natemp",sep="")))
			}
		}
	}

	eval(parse(text=paste("natemp <- ",depvar,sep="")))
	nay <- natemp[nahelp]
	eval(parse(text=paste(depvar," <- nay",sep="")))

	# end clean up missing values

	# defining useful functional forms -- mlog corrects a precision difficult - so log(x) does not return -Inf for x~0 

	mlog <- function(x) { ifelse(log(x) != "-Inf",log(x),-5000) }
	logit <- function(X) {(1 / (1 + exp(-X)))}

	"%d%" <- function(x,y) { (1 - (1-x)*(1-y)) };
	"%a%" <- function(x,y) { x*y } 
	
	# the probability structure is transformed into the likelihood function

	structure <- gsub("&","%a%",structure,extended=FALSE)
	structure <- gsub("|","%d%",structure,extended=FALSE)


	for (i in 1:length(vars)) {

		t <- rep(1,3)
		t[1] <- (1/as.integer(paste(regexpr(paste(vars[i],")",sep=""),structure,extended=FALSE))))
		t[2] <- (1/as.integer(paste(regexpr(paste(vars[i],"%",sep=""),structure,extended=FALSE))))
		t[3] <- (1/as.integer(paste(regexpr(paste(vars[i],"+",sep=""),structure,extended=FALSE))))
		t <- which.max(t)
		if (vars[i] == "cons") {

			if (method == "probit") {structure <- sub("(cons",paste("pnorm(b[",i,"]",sep=""),structure,extended=FALSE);t<-0}	
			if (method == "logit") {structure <- sub("(cons",paste("logit(b[",i,"]",sep=""),structure,extended=FALSE);t<-0}	
		}
	

		if (t == 1) {structure <- sub(paste(vars[i],")",sep=""),paste("(",vars[i],"*b[",i,"]))",sep=""),structure,extended=FALSE)}
		if (t == 2) {structure <- sub(paste(vars[i],"%",sep=""),paste("(",vars[i],"*b[",i,"])%",sep=""),structure,extended=FALSE)}
		if (t == 3) {structure <- sub(paste(vars[i],"+",sep=""),paste("(",vars[i],"*b[",i,"])+",sep=""),structure,extended=FALSE)}
		if (method == "probit") {structure <- sub(paste("((",vars[i],"*b[",i,"])",sep=""),paste("pnorm((",vars[i],"*b[",i,"])",sep=""),structure,extended=FALSE)}
        	if (method == "logit") {structure <- sub(paste("((",vars[i],"*b[",i,"])",sep=""),paste("logit((",vars[i],"*b[",i,"])",sep=""),structure,extended=FALSE)}
		    
       }




	eval(parse(text=paste("temp <-",depvar,sep="")))
	temp <- ((temp == 1)|(temp==0)|(temp=NA))
	if (mean(temp)!=1) {cat(paste("`depvar' must be binary \n"));break}
	eval(parse(text=paste(depvar,"<- as.integer(",depvar,")",sep="")))


	q<-paste("llik <- function(b) {sum(-1*(1-(",depvar,"))*mlog(1-",structure,")-",depvar,"*mlog(",structure,"),na.rm=TRUE)}",sep="")
	eval(parse(text=q,n=-1))

	# likelihood function is now saved as llik, defined as -1*llik


	# minimization of -likelihood using either optim, nlm, or genoud.  nlm maximization continues in a 'smart search'
	# for iteration time defined in `safety'

	tempfunc <- runif(length(vars))
	maxop <- ""
	if (maxoptions != "") {maxop <- paste(",",maxop,sep="")}
	if (optimizer=="nlm") {
		out <- list()
		out$minimum <- 100000000
		for (i in 1:safety) {
		temp <- paste("nlm(llik,tempfunc,hessian = TRUE, iterlim = 10000",maxop,")")
		tempout <- eval(parse(text=temp))
		if (tempout$minimum < out$minimum) {out <- tempout;tempfunc <- (-1*out$estimate)}
		if (tempout$minimum == out$minimum) {tempfunc <- runif(length(vars))}
		if (det(tempout$hessian) == 0) {tempfunc <- runif(length(vars))}
		if (tempout$minimum > out$minimum) {tempfunc <- (-1*tempout$estimate)}

		}
	}

	#note: optim output is converted to the same form as nlm output


	if (optimizer=="optim") {

	out <- optim(tempfunc,llik,hessian = TRUE, method="BFGS",control=list(maxit=500))
	out$estimate <- out$par
	out$gradient <- matrix()
	out$minimum <- out$value
	out$iterations <- out$counts[1]
	
	}
	
	if (optimizer=="genoud") {

	library("rgenoud")        
	temp <- paste("genoud(llik, nvars=",length(tempfunc),", BFGS=TRUE, hessian=TRUE, pop.size=",popsize,")", sep="")
	out <- eval(parse(text=temp))
	out$estimate <- out$par
	out$gradient <- out$gradients
        out$gradient <- matrix()
	out$minimum <- out$value
	out$iterations <- out$generations
	
	}

	# if bootstrap option is set to `FALSE', standard errors, etc. are derived from the hessian -
	# if hessian is noninvertible, all are set to 0

	if (bootstrap == FALSE) {

		if (det(out$hessian)!=0) {
			zs <- abs(out$estimate)/sqrt(diag(solve(out$hessian)))
			confl <- out$estimate-1.96*sqrt(diag(solve(out$hessian)))
			confh <- out$estimate+1.96*sqrt(diag(solve(out$hessian)))
			ses <- sqrt(diag(solve(out$hessian)))
		}

		if (det(out$hessian)==0) {
			zs <- 0
			confl <- 0
			confh <- 0
			ses <- 0
		}

	}

	
	# bootstrap routine proceeds by selecting indices with replacement for `bootsize' samples, and performing
	# likelihood maximization using a weighed average

	if (bootstrap == TRUE) {

		tempfunc <- out$estimate
		eval(parse(text=paste("size <- length(",depvar,")",sep="")))
		for (i in 1:length(vars)) {eval(parse(text=paste("vardraws",i," <- integer(bootsize)",sep="")))}
		for (i in 1:bootsize) {
			indices <- 1:size
			sampindices <- sample(indices, size, replace = TRUE, prob = NULL)
			weights <- integer()
			for (j in 1:size) {weights[j] <- sum(sampindices == j)}
			q<-paste("llik <- function(b) {ifelse((max(abs(b))<10),weighted.mean(-1*(1-(",depvar,"))*mlog(1-",structure,")-",depvar,"*mlog(",structure,"), weights, na.rm=FALSE),1000)}",sep="")	
			eval(parse(text=q,n=-1))
			temp <- paste("nlm(llik,tempfunc,hessian = FALSE, iterlim = 100",maxop,")")
			tempout <- eval(parse(text=temp))
			for (k in 1:length(vars)) {eval(parse(text=paste("vardraws",k,"[",i,"] <- tempout$estimate[",k,"]",sep="")))}
		}


		for (k in 1:length(vars)) {eval(parse(text=paste("vardraws",k,"<- sort(vardraws",k,")",sep="")))}

		confl <- integer(length(vars))
		confh <- integer(length(vars))
		for (k in 1:length(vars)) {eval(parse(text=paste("confl[",k,"] <- vardraws",k,"[ceiling(.05*bootsize)]",sep="")))}
		for (k in 1:length(vars)) {eval(parse(text=paste("confh[",k,"] <- vardraws",k,"[ceiling(.95*bootsize)]",sep="")))}
		ses <- (confh-confl)/4
		zs <- abs(out$estimate)/ses

	}


	# original data are restored in cases where missing values have been excluded

	eval(parse(text=paste(depvar," <- natemp",sep="")))
	modelvarsnum <- 0
	for (i in 1:length(vars)) {
		if (vars[i] != "cons") {
			if (sum(vars[i] == vars[1:i])==1){
				modelvarsnum <- modelvarsnum + 1
				j <- modelvarsnum
				eval(parse(text=paste("holds",j," -> natemp",sep="")))
				eval(parse(text=paste(vars[i]," <- natemp",sep="")))
			}
		}
	}

	eval(parse(text=paste(depvar," <- natemp",sep="")))


	# end restoration of original data

	final<-booltest(calc=structure, loglik=-1*out$minimum, vars=vars, coefs=out$estimate, ses=ses, iter=out$iterations,
hess=out$hessian, grad=out$gradient, 
zsc = zs, pz = 1-pnorm(abs(out$estimate),mean=0,sd=abs(ses)),
c95lo=confl, c95hi=confh, pstruc=pstructure,meth=method)
	return(final)  

} 


boolfirst <- function(object,gvar,instance = 0,range = 0) {

	# object must be output from boolean, instance is used only when a single variable is used more than once 
	# in the model, range can be used for user defined ranges.

	# redefining useful functional forms
	
	mlog <- function(x) { ifelse(log(x) != "-Inf",log(x),-5000) }
	logit <- function(X) {(1 / (1 + exp(-X)))}
	"%d%" <- function(x,y) { (1 - (1-x)*(1-y)) };
	"%a%" <- function(x,y) { x*y } 
	
	# import method (logit or probit) from previous run of boolean

	method = object@method
	vars <- object@Variables
	values <- object@Coefficients
	structure <- object@pstructure
	index <- integer()
	dpoints <- 50
	
	# making sure variable selected is actually in model

	for (i in 1:length(vars)) {index[i] <- match(gvar,vars[i],nomatch=0)}
	if (sum(index) ==0) {cat(paste("No such variable in model. \n"));break}
	
	# where is the variable of interest? then, define its range
	
	index <- match(1,index)
	if (mean(abs(range)) == 0) {

		temp <- paste("min(",vars[index],")",sep="")
		q123232 <- eval(parse(text=temp))
		temp <- paste("max(",vars[index],")",sep="")
		q223232 <- eval(parse(text=temp))
		range<-seq(q123232,q223232,by=(q223232-q123232)/dpoints)

	}	
	
	# set all other variables to mean values

	mvars <- numeric()
	cons <- 1
	for (i in 1:length(vars)) {

		temp <- paste("mvars[i] <- mean(",vars[i],")",sep="")
		eval(parse(text=temp))

	 }

	# deduce dependent variable, probability structure

	t <- as.integer(paste(regexpr("~",structure,extended=FALSE)))
	temp <- substr(structure,1,(t-1))
	depvar <- sub(" ","",temp,extended=FALSE)
	structure <- substr(structure,(t+1),nchar(structure))
	structure <- gsub("&","%a%",structure,extended=FALSE)
	structure <- gsub("|","%d%",structure,extended=FALSE)

	tempindex <- length(structure)
	holder <- "qghgh"
	for (i in 1:length(structure)) {holder <- sub("qghgh",paste(structure[i],"qghgh",sep=""),holder,extended=FALSE)}
	structure <- sub("qghgh","",holder,extended=FALSE)
	holder <- "qghgh"
	for (i in 1:length(depvar)) {holder <- sub("qghgh",paste(depvar[i],"qghgh",sep=""),holder,extended=FALSE)}
	depvar <- sub("qghgh","",holder,extended=FALSE)

	# convert probability structure to usuable functional form

	for (i in 1:length(vars)) {
	
		if (vars[i] != vars[index]) {
			t <- rep(1,3)
			t[1] <- (1/as.integer(paste(regexpr(paste(vars[i],")",sep=""),structure,extended=FALSE))))
			t[2] <- (1/as.integer(paste(regexpr(paste(vars[i],"%",sep=""),structure,extended=FALSE))))
			t[3] <- (1/as.integer(paste(regexpr(paste(vars[i],"+",sep=""),structure,extended=FALSE))))
			t <- which.max(t)	
			if (vars[i] == "cons") {

				if (method == "probit") {structure <- sub("(cons",paste("pnorm(",values[i],sep=""),structure,extended=FALSE);t<-0}	
				if (method == "logit") {structure <- sub("(cons",paste("logit(",values[i],sep=""),structure,extended=FALSE);t<-0}	
			}

			if (t == 1) {structure <- sub(paste(vars[i],")",sep=""),paste("(",mvars[i],"*",values[i],"))",sep=""),structure,extended=FALSE)}
			if (t == 2) {structure <- sub(paste(vars[i],"%",sep=""),paste("(",mvars[i],"*",values[i],")%",sep=""),structure,extended=FALSE)}
			if (t == 3) {structure <- sub(paste(vars[i],"+",sep=""),paste("(",mvars[i],"*",values[i],")+",sep=""),structure,extended=FALSE)}
			if (method == "probit") {structure <- sub(paste("((",mvars[i]*values[i],")",sep=""),paste("pnorm((",mvars[i]*values[i],")",sep=""),structure,extended=FALSE)}
        		if (method == "logit")  {structure <- sub(paste("((",mvars[i]*values[i],")",sep=""),paste("logit((",mvars[i]*values[i],")",sep=""),structure,extended=FALSE)}
	 
    		}


		if (vars[i] == vars[index]) {
			t[1] <- (1/as.integer(paste(regexpr(paste(vars[i],")",sep=""),structure,extended=FALSE))))
			t[2] <- (1/as.integer(paste(regexpr(paste(vars[i],"%",sep=""),structure,extended=FALSE))))
			t[3] <- (1/as.integer(paste(regexpr(paste(vars[i],"+",sep=""),structure,extended=FALSE))))
			t <- which.max(t)	
			if (t == 1) {structure <- sub(paste(vars[i],")",sep=""),paste("(gsxgsx*",values[i],"))",sep=""),structure,extended=FALSE)}
			if (t == 2) {structure <- sub(paste(vars[i],"%",sep=""),paste("(gsxgsx*",values[i],")%",sep=""),structure,extended=FALSE)}
			if (t == 3) {structure <- sub(paste(vars[i],"+",sep=""),paste("(gsxgsx*",values[i],")+",sep=""),structure,extended=FALSE)}
			if (method == "probit") {structure <- sub(paste("((gsxgsx*",values[i],")",sep=""),paste("pnorm((gsxgsx*",values[i],")",sep=""),structure,extended=FALSE)}
			if (method == "logit")  {structure <- sub(paste("((gsxgsx*",values[i],")",sep=""),paste("logit((gsxgsx*",values[i],")",sep=""),structure,extended=FALSE)}
		
 	        }
			
	}

	# write likelihood function, store to `nllik'	

	tempfunc <- paste("nllik <- function(gsxgsx) {(",structure,")}",sep="")
	eval(parse(text=tempfunc))

	# first difference values calculated, plotted

	llikpts <- integer()
	for (i in 1:length(range)) {llikpts[i] <- nllik(range[i])}
	plot(range,llikpts,type="l",xlab=paste("Value of",vars[index]),ylab="Probability of Event")	
	
	
}





boolprep <- function(form, depvar, a, b = "", c = "", d = "", e = "", f = "", g = "", h = "", i = "",j = "", k = "", l = "", m = "", n = "", o = "", p = "", q = "", r = "", s = "", t = "", u = "", v = "", w = "", x = "", y = "", z = "", constant = TRUE) {

	if (!is.character(form)) {form <- deparse(form)}
	form <- gsub(" ","",form,extended=FALSE)
	t <- as.integer(paste(regexpr("~",form,extended=FALSE)))
	form <- substr(form,(t+1),nchar(form))
	form <- paste("(",form,")",sep="")
	print(form)
	
	for (i in letters) {

		if (constant == TRUE) {form <- sub(paste("(",i,"&",sep=""),paste("((cons+",eval(parse(text=i)),")&",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("(",i,"|",sep=""),paste("((cons+",eval(parse(text=i)),")|",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("&",i,"&",sep=""),paste("&(cons+",eval(parse(text=i)),")&",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("&",i,"|",sep=""),paste("&(cons+",eval(parse(text=i)),")|",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("&",i,")",sep=""),paste("&(cons+",eval(parse(text=i)),"))",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("|",i,"&",sep=""),paste("|(cons+",eval(parse(text=i)),")&",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("|",i,"|",sep=""),paste("|(cons+",eval(parse(text=i)),")|",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("|",i,")",sep=""),paste("|(cons+",eval(parse(text=i)),"))",sep=""),form,extended=FALSE)	
					}
		
		if (constant == FALSE) {form <- sub(paste("(",i,"&",sep=""),paste("((cons+",eval(parse(text=i)),")&",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("(",i,"|",sep=""),paste("((",eval(parse(text=i)),")|",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("&",i,"&",sep=""),paste("&(",eval(parse(text=i)),")&",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("&",i,"|",sep=""),paste("&(",eval(parse(text=i)),")|",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("&",i,")",sep=""),paste("&(",eval(parse(text=i)),"))",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("|",i,"&",sep=""),paste("|(",eval(parse(text=i)),")&",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("|",i,"|",sep=""),paste("|(",eval(parse(text=i)),")|",sep=""),form,extended=FALSE)	
		 	form <- sub(paste("|",i,")",sep=""),paste("|(",eval(parse(text=i)),"))",sep=""),form,extended=FALSE)		
					}


			   }

	paste(depvar,"~(",form,")")

}




boolprof <- function(object,gvar,instance = 0,range = 0) {

	# define useful functional forms
	
	mlog <- function(x) { ifelse(log(x) != "-Inf",log(x),-5000) }
	logit <- function(X) {(1 / (1 + exp(-X)))}
	"%d%" <- function(x,y) { (1 - (1-x)*(1-y)) }
	"%a%" <- function(x,y) { x*y } 

	# import necessary items from `object' which must be of type `booltest'

	method = object@method
	vars <- object@Variables
	values <- object@Coefficients
	structure <- object@pstructure
	
	# looking for variable of interest, various checks to ensure variable of interest is in model
	# and that if the variable occurs more than once, user has specified which instance

	index <- integer()
	dpoints <- 50
	for (i in 1:length(vars)) {index[i] <- match(gvar,vars[i],nomatch=0)}
	if (sum(index) ==0) {cat(paste("No such variable in model. \n"));break}
	if ((sum(index) > 1) & (instance == 0) ) {cat(paste("Variable occurs more than once in model. \n"));
			     cat(paste("Specify `instance' in function call. \n"));break}
	if (sum(index) == 1) {index <- match(1,index)}
	if ((sum(index) > 1) & (instance != 0) ) {cindex <- cumsum(index);
						  index <- match(instance,cindex,nomatch=0)}
	if (index==0) {cat(paste("The variable does not occur so many times. \n Reduce `instance in function call. \n"));break}
	if (sum(abs(range)) == 0) {range <- seq(values[index]-object@StandardErrors[index],values[index]+object@StandardErrors[index],by=(2*object@StandardErrors[index])/dpoints)}
	
	# deduce dependent variable, probability structure from items in booltest	

	t <- as.integer(paste(regexpr("~",structure,extended=FALSE)))
	temp <- substr(structure,1,(t-1))
	depvar <- sub(" ","",temp,extended=FALSE)	
	structure <- substr(structure,(t+1),nchar(structure))
	structure <- gsub("&","%a%",structure,extended=FALSE)
	structure <- gsub("|","%d%",structure,extended=FALSE)

	tempindex <- length(structure)
	holder <- "qghgh"
	for (i in 1:length(structure)) {holder <- sub("qghgh",paste(structure[i],"qghgh",sep=""),holder,extended=FALSE)}
	structure <- sub("qghgh","",holder,extended=FALSE)
	holder <- "qghgh"
	for (i in 1:length(depvar)) {holder <- sub("qghgh",paste(depvar[i],"qghgh",sep=""),holder,extended=FALSE)}
	depvar <- sub("qghgh","",holder,extended=FALSE)


	# transform probabilty structure into usuable functional form

	for (i in 1:length(vars)) {
	
		if (i != index) {
			t <- rep(1,3)
			t[1] <- (1/as.integer(paste(regexpr(paste(vars[i],")",sep=""),structure,extended=FALSE))))
			t[2] <- (1/as.integer(paste(regexpr(paste(vars[i],"%",sep=""),structure,extended=FALSE))))
			t[3] <- (1/as.integer(paste(regexpr(paste(vars[i],"+",sep=""),structure,extended=FALSE))))
			t <- which.max(t)	
			if (vars[i] == "cons") {
				if (method == "probit") {structure <- sub("(cons",paste("pnorm(",values[i],sep=""),structure,extended=FALSE);t<-0}	
				if (method == "logit") {structure <- sub("(cons",paste("logit(",values[i],sep=""),structure,extended=FALSE);t<-0}	
			}
			if (t == 1) {structure <- sub(paste(vars[i],")",sep=""),paste("(",vars[i],"*",values[i],"))",sep=""),structure,extended=FALSE)}
			if (t == 2) {structure <- sub(paste(vars[i],"%",sep=""),paste("(",vars[i],"*",values[i],")%",sep=""),structure,extended=FALSE)}
			if (t == 3) {structure <- sub(paste(vars[i],"+",sep=""),paste("(",vars[i],"*",values[i],")+",sep=""),structure,extended=FALSE)}
			if (method == "probit") {structure <- sub(paste("((",vars[i],"*",values[i],")",sep=""),paste("pnorm((",vars[i],"*",values[i],")",sep=""),structure,extended=FALSE)}
   			if (method == "logit")  {structure <- sub(paste("((",vars[i],"*",values[i],")",sep=""),paste("logit((",vars[i],"*",values[i],")",sep=""),structure,extended=FALSE)}
	   	}

	
		if (i == index) {
			t[1] <- (1/as.integer(paste(regexpr(paste(vars[i],")",sep=""),structure,extended=FALSE))))
			t[2] <- (1/as.integer(paste(regexpr(paste(vars[i],"%",sep=""),structure,extended=FALSE))))
			t[3] <- (1/as.integer(paste(regexpr(paste(vars[i],"+",sep=""),structure,extended=FALSE))))
			t <- which.max(t)	
			if (vars[i] == "cons") {
				if (method == "probit") {structure <- sub("(cons","pnorm(gsxgsx",structure,extended=FALSE);t<-0}	
				if (method == "logit") {structure <- sub("(cons","logit(gsxgsx",structure,extended=FALSE);t<-0}	
			}
		if (t == 1) {structure <- sub(paste(vars[i],")",sep=""),paste("(",vars[i],"*gsxgsx))",sep=""),structure,extended=FALSE)}
		if (t == 2) {structure <- sub(paste(vars[i],"%",sep=""),paste("(",vars[i],"*gsxgsx)%",sep=""),structure,extended=FALSE)}
		if (t == 3) {structure <- sub(paste(vars[i],"+",sep=""),paste("(",vars[i],"*gsxgsx)+",sep=""),structure,extended=FALSE)}
		if (method == "probit") {structure <- sub(paste("((",vars[i],"*gsxgsx)",sep=""),paste("probit((",vars[i],"*gsxgsx)",sep=""),structure,extended=FALSE)}
		if (method == "logit") {structure <- sub(paste("((",vars[i],"*gsxgsx)",sep=""),paste("logit((",vars[i],"*gsxgsx)",sep=""),structure,extended=FALSE)}	

		}
	}

	# write likelihood function
	
	tempfunc <- paste("nllik <- function(gsxgsx) {sum((1-(",depvar,"))*mlog(1-",structure,")+",depvar,"*mlog(",structure,"))}",sep="")
	eval(parse(text=tempfunc))

	# evalute likelihood function over two standard deviations around the mean, plot these evaluations

	llikpts <- integer()
	for (i in 1:length(range)) {llikpts[i] <- nllik(range[i])}
	plot(range,llikpts,type="l",xlab=paste("Coefficient of",vars[index]),ylab="Value of Likelihood Function")	
		
}



