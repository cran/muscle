#include "muscle.h"
#include <stdio.h>

#include <R.h>

static char szOnExceptionMessage[] =
	{
	"\nFatal error, exception caught.\n"
	};

void OnException()
	{
	Rprintf("%s", szOnExceptionMessage);
	Log("%s", szOnExceptionMessage);
	Log("Finished %s\n", GetTimeAsStr());
	error("error in muscle\n");
	}


