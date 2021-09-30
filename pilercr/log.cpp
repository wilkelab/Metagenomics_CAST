#include "pilercr.h"

static FILE *g_fLog = 0;

void SetLog()
	{
	bool Append = false;
	const char *FileName = ValueOpt("log");
	if (0 == FileName)
		{
		FileName = ValueOpt("loga");
		if (0 == FileName)
			return;
		Append = true;
		}

	g_fLog = fopen(FileName, Append ? "a" : "w");
	if (NULL == g_fLog)
		{
		fprintf(stderr, "\n*** FATAL ERROR *** Cannot open log file\n");
		perror(FileName);
		exit(1);
		}
	setbuf(g_fLog, 0);
	}

void Log(const char Format[], ...)
	{
	if (0 == g_fLog)
		return;

	char Str[4096];
	va_list ArgList;
	va_start(ArgList, Format);
	vsprintf(Str, Format, ArgList);
	fprintf(g_fLog, "%s", Str);
	if (g_FlushLog)
		fflush(g_fLog);
	}

void Out(const char Format[], ...)
	{
	if (0 == g_fOut)
		return;

	char Str[4096];
	va_list ArgList;
	va_start(ArgList, Format);
	vsprintf(Str, Format, ArgList);
	fprintf(g_fOut, "%s", Str);
	Log("%s", Str);
	}
