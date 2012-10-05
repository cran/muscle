//@@TODO reconcile /muscle with /muscle3.6

// Author: Robert C. Edgar
// Ported into R by Alex T. Kalinka (alex.t.kalinka@gmail.com)


#include "muscle.h"
#include <stdio.h>
#ifdef	WIN32
#include <windows.h>	// for SetPriorityClass()
#include <io.h>			// for isatty()
#else
#include <unistd.h>		// for isatty()
#endif


#include <R.h>


const char *MUSCLE_LONG_VERSION	= "MUSCLE v" SHORT_VERSION "."
#include "svnversion.h"
" by Robert C. Edgar";


int g_argc;
char **g_argv;



extern "C" {



void muscleR(int *argc, char **argv)
	{

	int nargs = *argc;

#if	WIN32
// Multi-tasking does not work well in CPU-bound
// console apps running under Win32.
// Reducing the process priority allows GUI apps
// to run responsively in parallel.
	SetPriorityClass(GetCurrentProcess(), BELOW_NORMAL_PRIORITY_CLASS);
#endif
	g_argc = nargs;
	g_argv = argv;

	SetNewHandler();
	SetStartTime();
	ProcessArgVect(nargs, argv);
	SetParams();
	SetLogFile();

	//extern void TestSubFams(const char *);
	//TestSubFams(g_pstrInFileName);
	//return 0;

	if (g_bVersion)
		{
		Rprintf("%s\n", MUSCLE_LONG_VERSION);
		return;
		}

	if (!g_bQuiet)
		Credits();

	if (MissingCommand() && isatty(0))
		{
		Usage();
		return;
		}

	if (g_bCatchExceptions)
		{
		try
			{
			Run();
			}
		catch (...)
			{
			OnException();
			return;
			}
		}
	else
		Run();

	return;
	}


}

