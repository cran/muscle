#include "muscle.h"
#include "msa.h"
#include "textfile.h"
#include <time.h>

#include <R.h>

MSA *ptrBestMSA;
static const char *pstrOutputFileName;

void SetOutputFileName(const char *out)
	{
	pstrOutputFileName = out;
	}

void SetCurrentAlignment(MSA &msa)
	{
	ptrBestMSA = &msa;
	}

void SaveCurrentAlignment()
	{
	static bool bCalled = false;
	if (bCalled)
		{
		Rprintf(
		  "\nRecursive call to SaveCurrentAlignment, giving up attempt to save.\n");
		return;
		}

	if (0 == ptrBestMSA)
		{
		Rprintf("\nAlignment not completed, cannot save.\n");
		Log("Alignment not completed, cannot save.\n");
		return;
		}

	if (0 == pstrOutputFileName)
		{
		Rprintf("\nOutput file name not specified, cannot save.\n");
		return;
		}

	Rprintf("\nSaving current alignment ...\n");

	TextFile fileOut(pstrOutputFileName, true);
	ptrBestMSA->ToFASTAFile(fileOut);

	Rprintf("Current alignment saved to \"%s\".\n", pstrOutputFileName);
	Log("Current alignment saved to \"%s\".\n", pstrOutputFileName);
	}

void CheckMaxTime()
	{
	if (0 == g_ulMaxSecs)
		return;

	time_t Now = time(0);
	time_t ElapsedSecs = Now - GetStartTime();
	if (ElapsedSecs <= (time_t) g_ulMaxSecs)
		return;

	Log("Max time %s exceeded, elapsed seconds = %ul\n",
	  MaxSecsToStr(), ElapsedSecs);

	SaveCurrentAlignment();
	return;
	}
