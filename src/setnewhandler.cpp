#include "muscle.h"
#include <stdio.h>
#include <new>

#include <R.h>

const int ONE_MB = 1024*1024;
const size_t RESERVE_BYTES = 8*ONE_MB;
static void *EmergencyReserve = 0;

void OnOutOfMemory()
	{
	free(EmergencyReserve);
	Rprintf("\n*** OUT OF MEMORY ***\n");
	Rprintf("Memory allocated so far %g MB\n", GetMemUseMB());
	extern MSA *ptrBestMSA;
	if (ptrBestMSA == 0)
		Rprintf("No alignment generated\n");
	else
		SaveCurrentAlignment();
	return;
	}

void SetNewHandler()
	{
	EmergencyReserve = malloc(RESERVE_BYTES);
	std::set_new_handler(OnOutOfMemory);
	}
