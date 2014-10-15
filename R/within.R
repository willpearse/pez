##' @export
within.comparative.comm <- function(data, expr, ...){
	
	parent <- parent.frame()
	e <- evalq(environment(), data, parent)
	eval(substitute(expr), e)
	
	cc <- as.list(e)
	class(cc) <- class(data)
        return(cc)
}


## FIXME: is this dangerous somehow?

