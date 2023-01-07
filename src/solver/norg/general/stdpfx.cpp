/*
code developed by
    Rong-Qiang He (rqhe@ruc.edu.cn, RUC, China)
date 2013 - 2017
*/

#include "stdpfx.h"

#ifdef _MSC_VER
// date and time
std::string present()
{
	struct tm newtime;
	time_t seconds;

	time(&seconds);
	localtime_s(&newtime, &seconds);

	char timebuf[26];
	sprintf_s(timebuf, 26, "%d-%02d-%02d %02d:%02d:%02d",
		newtime.tm_year + 1900,
		newtime.tm_mon + 1,
		newtime.tm_mday,
		newtime.tm_hour,
		newtime.tm_min,
		newtime.tm_sec);
	return std::string(timebuf);
}
// date and time with default format
std::string present_time_default_format()
{
	struct tm newtime;
	time_t seconds;

	time(&seconds);
	localtime_s(&newtime, &seconds);

	char timebuf[26];
	asctime_s(timebuf, 26, &newtime);
	return std::string(timebuf);
}
#else
// date and time
std::string present()
{
	struct tm newtime;
	time_t seconds;

	time(&seconds);
	localtime_r(&seconds, &newtime);

	char timebuf[26];
	snprintf(timebuf, 26, "%d-%02d-%02d %02d:%02d:%02d",
		newtime.tm_year + 1900,
		newtime.tm_mon + 1,
		newtime.tm_mday,
		newtime.tm_hour,
		newtime.tm_min,
		newtime.tm_sec);
	return std::string(timebuf);
}
// date and time with default format
std::string present_time_default_format()
{
	struct tm newtime;
	time_t seconds;

	time(&seconds);
	localtime_r(&seconds, &newtime);

	char timebuf[26];
	asctime_r(&newtime, timebuf);
	return std::string(timebuf);
}
#endif
