/*
 * AudSWIPEP.h
 *
 *  Created on: Apr 24, 2013
 *      Author: saul
 */

#ifndef AUDSWIPEP_H_
#define AUDSWIPEP_H_
#include "includes.h"
#include "Clocks.h"
#include "MPIMatlabCommunicator.h"
#include "SPUtilities.h"
#include "AuditiveSystem.h"
struct timespec start, end, startTotal, endTotal;
clocksArray clocks, cycleClocks, totalTime;
int myid;
double soundLength;
#endif
