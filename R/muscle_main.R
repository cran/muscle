#####################################################################################
# R wrapper function(s) for Edgar's multiple sequence alignment (MUSCLE) software.  #
# see:                                                                              #
# Edgar, R.C. (2004)                                                                #
#	MUSCLE: multiple sequence alignment with high accuracy and high throughput. #
#	Nucleic Acids Res 32, 1792-1797.                                            #
# and:                                                                              #
# http://www.drive5.com/muscle/muscle_userguide3.8.html                             #
#                                                                                   #
# Original author of MUSCLE algorithm: Robert C. Edgar                              #
# Ported into R by Alex T. Kalinka (alex.t.kalinka@gmail.com)                       #
#                                                                                   #
#####################################################################################



print.muscle <- function(x, ...)
	# S3 method for generic function "print".
	# x is a "muscle" object.
	{
	args <- list(...)
	if(length(args)>0){
		mf <- match("from",names(args))
		mt <- match("to",names(args))
		mn <- match("num",names(args))
		ms <- match("seqs",names(args))
		if(!is.na(mf)){
			if(is.na(mt)){
				from <- args$from
				to <- from + 40
				if(from > x$length){
					from <- x$length-40
					}
				if(to > x$length){
					to <- x$length
					}
			}else{
				from <- args$from
				to <- args$to
				if(from > x$length){
					from <- x$length-40
					}
				if(to > x$length){
					to <- x$length
					}
				}
			}
		if(!is.na(mn)){
			num <- 1:args$num
		}else{
			num <- 1:nrow(x$seqs)
			}
		if(!is.na(ms)){
			num <- match(args$seqs,x$seqs[,1])
			if(length(which(is.na(num)==TRUE))>0){
				num <- 1:nrow(x$seqs)
				}
			}
	}else{
		from <- 1
		to <- 40
		num <- 1:nrow(x$seqs)
		}
	summ <- NULL
	cmx <- max(sapply(x$seqs[num,1],nchar))
	wd <- getOption("width")-5-cmx-4
	tot <- to-from+1
	if(tot > wd){
		il <- ceiling((tot/wd))
		len <- wd
		lt <- len+from-1; plt <- lt
	}else{
		il <- 1
		len <- to-from+1
		lt <- to; plt <- lt
		}
	lf <- from
	for(i in 1:il){ # interleave sequences.
		k <- 1
		for(j in num){
			add <- cmx - nchar(x$seqs[j,1])
			if(k==1){
				if((lt-lf+1) < len){
					jl <- lt-lf+1
				}else{
					jl <- len
					}
				nc <- nchar(as.character(lf))
				summ <- append(summ, c(rep(" ",(cmx+4-nc)),lf,rep(" ",jl),lt,"\n"))
				}
			summ <- append(summ, c(rep(" ",add), as.character(x$seqs[j,1]), "    ", as.character(substr(x$seqs[j,2],lf,lt)),"\n"))
			k <- k+1
			}
		lf <- lt+1
		lt <- lt+len
		if(lt > to){
			lt <- to
			}
		summ <- append(summ,"\n")
		}
	summ <- paste(summ, collapse="")
	cat("\n# ",x$num," sequences with ",x$length," positions\n\n","# Position ",from," to ",to,":\n\n",summ,"\n",sep="")
	}


.expand.tilde <- function(file)
	{
	uu <- unlist(strsplit(file,""))
	gg <- grep("~",uu)
	if(length(gg) == 1){ # Need to expand tilde?
		if(capabilities("cledit")){ # readline available?
			if(gg == 1){
				file <- path.expand(file)
			}else{
				stop("the tilde character (~) must be the first character in the file path\n")
				}
		}else{
			stop("readline not available: absolute file path required\n i.e. do not use the tilde character (~)\n")
			}
	}else if(length(gg) > 1){
		stop("there must be only one tilde character (~) in the file path\n")
		}
	return(file)
	}


muscle <- function(seqs, out = NULL, quiet = FALSE, ...)
	{
	# seqs is either a character vector specifying a fasta file of sequences or an object of class "fasta".
	# options for alternative output are listed at the drive5 HTML help page.

	args <- list(...)
	arg.v <- NULL
	write <- FALSE

	if(is.character(seqs)){
		seqs <- .expand.tilde(seqs)
		if(file.access(seqs) == -1){
			stop(cat("\nfile not found: \"",seqs,"\"\n",sep=""))
			}
		arg.v <- append(arg.v,c("in",as.character(seqs)))
	}else if(class(seqs)=="fasta" || class(seqs)=="muscle"){
		write.fasta(seqs, file="temp.fa")
		arg.v <- append(arg.v,c("in",as.character("temp.fa")))
		write <- TRUE
	}else{
		stop("input \'seqs\' should be a file name, or an object of class \'muscle\' or \'fasta\'\n")
		}
	if(is.null(out)){
		arg.v <- append(arg.v,c("out","temp.afa"))
	}else{
		arg.v <- append(arg.v,c("out",as.character(out)))
		}
	if(quiet){
		arg.v <- append(arg.v,as.character("quiet"))
		}

	if(length(args) > 0){
		for(i in 1:length(args)){
			if(is.logical(args[[i]])){
				# Is a flag.
				arg.v <- append(arg.v, as.character(names(args)[i]))
			}else{ 
				# Is an option.
				arg.v <- append(arg.v, c(as.character(names(args)[i]),as.character(args[[i]])))
				}
			}
		}
	nargs <- as.integer(length(arg.v))

	stuff <- .C("muscleR", nargs, as.character(arg.v))
	if(write){
		file.remove("temp.fa")
		}

	if(is.null(out)){
		.C("read_fasta", file = as.character("temp.afa"))
		aln <- read.table("temp.rafa", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
		file.remove("temp.afa","temp.rafa")
		ret <- list()
		ret$seqs <- aln
		ret$length <- nchar(as.character(aln[1,2]))
		ret$num <- nrow(aln)
		class(ret) <- "muscle"
		return(ret)
	}else{
		return(invisible())
		}

	}


read.fasta <- function(file)
	{
	file <- .expand.tilde(file)
	.C("read_fasta", file = as.character(file))
	aln <- read.table("temp.rafa", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
	file.remove("temp.rafa")
	ret <- list()
	ret$seqs <- aln
	ret$num <- nrow(aln)
	class(ret) <- "fasta"
	return(ret)
	}


write.fasta <- function(aln, file)
	{
	# aln is an object of class "muscle" or "fasta".
	seqs <- NULL
	if(class(aln)=="muscle" || class(aln)=="fasta"){
		num <- nrow(aln$seqs)*2
		#len <- nchar(aln$seqs[1,2])
		for(i in 1:(num/2)){
			seqs <- append(seqs, c(as.character(aln$seqs[i,1]), as.character(aln$seqs[i,2])))
			}
	}else{
		stop("input must be an object of class \'muscle\' or \'fasta\'\n")
		}
	file <- .expand.tilde(file)
	out <- .C("write_fasta", as.character(seqs), as.character(file), as.integer(num))
	return(invisible())
	}



